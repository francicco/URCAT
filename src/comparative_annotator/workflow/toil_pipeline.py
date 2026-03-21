from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from toil.common import Toil
from toil.job import Job

from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.missing.consensus import (
    build_consensus_missing_transcript,
    choose_missing_locus_strand,
    cluster_projected_transcripts,
)
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus
from comparative_annotator.projection.reconstruct import reconstruct_projected_transcripts
from comparative_annotator.workflow.progressive import (
    compute_reference_order,
    extract_missing_locus_payloads,
    locus_overlaps_any_interval,
    pick_next_reference_from_order,
)

from comparative_annotator.workflow.orthology_edges import build_target_edge_evidence

from comparative_annotator.workflow.projection_table import write_projection_evidence_table
from comparative_annotator.workflow.novel_annotation_table import write_novel_annotations_table
from comparative_annotator.workflow.round_metrics import write_round_metrics_table
from comparative_annotator.workflow.workflow_paths import get_round_dir

from comparative_annotator.workflow.analysis_tables import (
    write_all_analysis_tables_for_round,
    write_all_analysis_tables_for_target,
)

from comparative_annotator.workflow.final_gff3 import write_all_final_species_gff3

@dataclass
class FrontierSeedTranscript:
    transcript_id: str
    species: str
    seqid: str
    strand: str
    exons: list


def read_json(path):
    path = Path(path)
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def chunked(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def get_species_list(species_csv: str):
    return [x.strip() for x in species_csv.split(",") if x.strip()]


def append_unique_preserve_order(items, value):
    out = list(items)
    if value not in out:
        out.append(value)
    return out


def load_transcripts_for_species(annotation_dir: str, annotation_suffix: str, species: str):
    gff_path = (Path(annotation_dir) / f"{species}{annotation_suffix}").resolve()
    return load_gff3(str(gff_path), species=species)


def load_all_transcripts(annotation_dir: str, annotation_suffix: str, species_list: list[str]):
    return {
        sp: load_transcripts_for_species(annotation_dir, annotation_suffix, sp)
        for sp in species_list
    }


def build_all_species_loci(transcripts_by_species):
    species_loci = {}
    for sp, txdict in transcripts_by_species.items():
        species_loci[sp] = build_species_loci(list(txdict.values()), species=sp)
    return species_loci


def get_hal_tree_newick(hal_path: str) -> str:
    result = subprocess.run(
        ["halStats", hal_path],
        capture_output=True,
        text=True,
        check=True,
    )
    for line in result.stdout.splitlines():
        line = line.strip()
        if line.endswith(";") and "(" in line:
            return line
    raise RuntimeError("Could not extract species tree from halStats output")


def build_consensus_seed(item):
    return FrontierSeedTranscript(
        transcript_id=(
            f"URCAT_CONSENSUS:{item['species']}:{item['seqid']}:"
            f"{item['start']}-{item['end']}"
        ),
        species=item["species"],
        seqid=item["seqid"],
        strand=item["strand"],
        exons=item["exons"],
    )


def collect_explained_locus_ids_from_round(round_merged_json):
    explained = {}

    for target_block in round_merged_json["targets"]:
        target_species = target_block["target_species"]
        explained_ids = set()

        for r in target_block["results"]:
            if r.get("status") != "ok":
                continue

            for sp, locus_id in r.get("primary", {}).items():
                if sp == target_species and locus_id:
                    explained_ids.add(locus_id)

            for sp, locus_ids in r.get("alternatives", {}).items():
                if sp == target_species:
                    explained_ids.update(locus_ids)

        explained[target_species] = explained_ids

    return explained


def collect_projected_transcript_spans_for_species(
    transcripts_by_species,
    source_species_list,
    target_species,
    hal,
):
    spans = []

    for source_species in source_species_list:
        if source_species == target_species:
            continue

        for _, seed in transcripts_by_species[source_species].items():
            projected_exon_blocks = []

            for exon_start, exon_end in seed.exons:
                intervals = hal.project_interval(
                    source_species=seed.species,
                    target_species=target_species,
                    seqid=seed.seqid,
                    start=exon_start,
                    end=exon_end,
                    strand=seed.strand,
                    source_transcript=seed.transcript_id,
                )
                projected_exon_blocks.append(intervals)

            pts = reconstruct_projected_transcripts(seed, projected_exon_blocks)
            for pt in pts:
                spans.append(
                    {
                        "seqid": pt.seqid,
                        "start": pt.start,
                        "end": pt.end,
                        "strand": pt.strand,
                        "source_species": source_species,
                        "source_transcript": pt.source_transcript,
                    }
                )

    return spans


def build_round_summary(round_merged, decision):
    reference_species = round_merged["reference_species"]
    round_id = round_merged["round_id"]

    targets_summary = {}

    for target_block in round_merged["targets"]:
        target_species = target_block["target_species"]
        results = target_block["results"]

        n_processed = len(results)
        n_ok = sum(1 for r in results if r.get("status") == "ok")
        n_error = sum(1 for r in results if r.get("status") == "error")

        n_primary = sum(
            1
            for r in results
            if r.get("status") == "ok" and bool(r.get("primary"))
        )

        n_alternative_only = sum(
            1
            for r in results
            if r.get("status") == "ok"
            and not bool(r.get("primary"))
            and bool(r.get("alternatives"))
        )

        n_missing = sum(
            1
            for r in results
            if r.get("status") == "ok" and bool(r.get("missing_annotations"))
        )

        n_strand_conflict = sum(
            1
            for r in results
            if r.get("status") == "ok" and bool(r.get("strand_conflicts"))
        )

        n_new_consensus = len(decision.get("new_consensus_by_species", {}).get(target_species, []))
        n_orphan_loci = len(decision.get("orphan_loci_by_species", {}).get(target_species, []))
        n_pending_frontier = len(decision.get("pending_frontiers_by_species", {}).get(target_species, []))

        targets_summary[target_species] = {
            "n_processed": n_processed,
            "n_ok": n_ok,
            "n_error": n_error,
            "n_primary": n_primary,
            "n_alternative_only": n_alternative_only,
            "n_missing": n_missing,
            "n_strand_conflict": n_strand_conflict,
            "n_new_consensus": n_new_consensus,
            "n_orphan_loci": n_orphan_loci,
            "n_pending_frontier": n_pending_frontier,
        }

    return {
        "round_id": round_id,
        "seed_species": decision["seed_species"],
        "reference_species": reference_species,
        "reference_order": decision["reference_order"],
        "used_reference_species": decision["used_reference_species"],
        "next_reference_species": decision["next_reference_species"],
        "stop": decision["stop"],
        "targets": targets_summary,
    }


def finalize_round_outputs(output_dir: str, round_id: int) -> None:
    round_dir = get_round_dir(output_dir, round_id)
    write_all_analysis_tables_for_round(round_dir)

def write_round_summary(job, workdir, round_merged_path, decision_path):
    round_merged = read_json(round_merged_path)
    decision = read_json(decision_path)

    summary = build_round_summary(round_merged, decision)

    round_id = summary["round_id"]
    reference_species = summary["reference_species"]

    ref_dir = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
    )
    round_dir = Path(workdir) / "rounds" / f"round_{round_id:03d}"

    json_paths = [ref_dir / "summary.json", round_dir / "summary.json"]
    tsv_paths = [ref_dir / "summary.tsv", round_dir / "summary.tsv"]

    for json_path in json_paths:
        write_json(json_path, summary)

    header = [
        "round_id",
        "seed_species",
        "reference_species",
        "target_species",
        "n_processed",
        "n_ok",
        "n_error",
        "n_primary",
        "n_alternative_only",
        "n_missing",
        "n_strand_conflict",
        "n_new_consensus",
        "n_orphan_loci",
        "n_pending_frontier",
        "next_reference_species",
        "stop",
    ]

    for tsv_path in tsv_paths:
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(tsv_path, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for target_species, stats in summary["targets"].items():
                fh.write(
                    "	".join(
                        [
                            str(summary["round_id"]),
                            str(summary["seed_species"]),
                            str(summary["reference_species"]),
                            str(target_species),
                            str(stats["n_processed"]),
                            str(stats["n_ok"]),
                            str(stats["n_error"]),
                            str(stats["n_primary"]),
                            str(stats["n_alternative_only"]),
                            str(stats["n_missing"]),
                            str(stats["n_strand_conflict"]),
                            str(stats["n_new_consensus"]),
                            str(stats["n_orphan_loci"]),
                            str(stats["n_pending_frontier"]),
                            str(summary["next_reference_species"]),
                            str(summary["stop"]),
                        ]
                    )
                    + "\n"
                )

    return str(json_paths[0])


def write_seed_frontier(job, workdir, annotation_dir, annotation_suffix, seed_species):
    tx = load_transcripts_for_species(annotation_dir, annotation_suffix, seed_species)
    out = {
        "round_id": 0,
        "reference_species": seed_species,
        "frontier_kind": "seed_all_transcripts",
        "transcript_ids": sorted(tx.keys()),
    }

    path = (
        Path(workdir)
        / "rounds"
        / "round_000"
        / f"ref_{seed_species}"
        / "frontier.json"
    )
    write_json(path, out)
    return str(path)


def write_species_frontier_from_pending(
    job,
    workdir,
    pending_frontiers_by_species,
    reference_species,
    round_id,
):
    pending = pending_frontiers_by_species.get(reference_species, [])

    frontier = {
        "round_id": round_id,
        "reference_species": reference_species,
        "frontier_kind": "mixed_pending",
        "items": pending,
    }

    path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / "frontier.json"
    )
    write_json(path, frontier)
    return str(path)


def write_manifest(
    job,
    workdir,
    frontier_path,
    reference_species,
    target_species_list,
    batch_size,
):
    frontier = read_json(frontier_path)
    round_id = frontier["round_id"]

    if frontier["frontier_kind"] == "seed_all_transcripts":
        items = [
            {"kind": "native_transcript", "transcript_id": tid}
            for tid in frontier["transcript_ids"]
        ]
    elif frontier["frontier_kind"] == "mixed_pending":
        items = frontier["items"]
    else:
        raise ValueError(f"Unknown frontier kind: {frontier['frontier_kind']}")

    batches = list(chunked(items, batch_size))

    jobs = []
    for target in target_species_list:
        for batch_id, batch_items in enumerate(batches):
            batch_file = (
                Path(workdir)
                / "rounds"
                / f"round_{round_id:03d}"
                / f"ref_{reference_species}"
                / f"target_{target}"
                / f"batch_{batch_id:03d}.items.json"
            )
            write_json(batch_file, {"items": batch_items})

            jobs.append(
                {
                    "round_id": round_id,
                    "target_species": target,
                    "batch_id": batch_id,
                    "batch_file": str(batch_file),
                    "n_items": len(batch_items),
                }
            )

    manifest = {
        "round_id": round_id,
        "reference_species": reference_species,
        "batch_size": batch_size,
        "jobs": jobs,
    }

    manifest_path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / "manifest.json"
    )
    write_json(manifest_path, manifest)
    return str(manifest_path)


def run_project_batch(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    round_id,
    reference_species,
    target_species,
    batch_items_path,
    batch_id,
):
    species_list = get_species_list(species_csv)
    transcripts_by_species = load_all_transcripts(annotation_dir, annotation_suffix, species_list)
    species_loci = build_all_species_loci(transcripts_by_species)
    hal = HALAdapter(str(Path(hal_path).resolve()))

    batch_info = read_json(batch_items_path)
    items = batch_info["items"]

    results = []

    for item in items:
        try:
            if item["kind"] == "native_transcript":
                source_label = item["transcript_id"]
                seed = transcripts_by_species[reference_species][item["transcript_id"]]

            elif item["kind"] == "urcat_consensus":
                source_label = ",".join(item.get("source_transcripts", [])) or "URCAT_CONSENSUS"
                seed = build_consensus_seed(item)

            elif item["kind"] == "orphan_native_locus":
                tx_id = item["transcripts"][0]
                source_label = tx_id
                seed = transcripts_by_species[reference_species][tx_id]

            else:
                raise ValueError(f"Unknown seed kind: {item['kind']}")

            clocus = infer_comparative_locus(
                seed_transcript=seed,
                target_species=target_species,
                hal_adapter=hal,
                species_loci=species_loci,
                transcripts_by_species=transcripts_by_species,
            )

            result = {
                "source_species": reference_species,
                "source_transcript": source_label,
                "target_species": target_species,
                "seed_kind": item["kind"],
                "status": "ok",
                "primary": clocus.primary,
                "alternatives": clocus.alternatives,
                "missing_annotations": clocus.missing_annotations,
                "strand_conflicts": clocus.strand_conflicts,
                "primary_transcripts": getattr(clocus, "primary_transcripts", {}),
                "alternative_transcripts": getattr(clocus, "alternative_transcripts", {}),
            }

        except Exception as e:
            result = {
                "source_species": reference_species,
                "source_transcript": item.get("transcript_id", item.get("locus_id", "unknown")),
                "target_species": target_species,
                "seed_kind": item.get("kind", "unknown"),
                "status": "error",
                "error": str(e),
            }

        results.append(result)

    out = {
        "round_id": round_id,
        "reference_species": reference_species,
        "target_species": target_species,
        "batch_id": batch_id,
        "n_items": len(items),
        "results": results,
    }

    out_path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / f"target_{target_species}"
        / f"batch_{batch_id:03d}.json"
    )
    write_json(out_path, out)
    return str(out_path)


def merge_target_results(job, workdir, round_id, reference_species, target_species, batch_result_paths):
    batch_results = [read_json(p) for p in batch_result_paths]

    merged_results = []
    for batch in batch_results:
        merged_results.extend(batch["results"])

    merged = {
        "round_id": round_id,
        "reference_species": reference_species,
        "target_species": target_species,
        "n_batches": len(batch_results),
        "n_results": len(merged_results),
        "results": merged_results,
    }

    out_path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / f"target_{target_species}"
        / "merged.json"
    )
    write_json(out_path, merged)
    return str(out_path)

def run_target_edge_evidence(
    job,
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    hal_path: str,
    species_csv: str,
    merged_target_path: str,
) -> str:
    edge_json_path = build_target_edge_evidence(
        workdir=workdir,
        annotation_dir=annotation_dir,
        annotation_suffix=annotation_suffix,
        hal_path=hal_path,
        species_csv=species_csv,
        merged_target_path=merged_target_path,
    )
    target_dir = Path(edge_json_path).parent
    write_all_analysis_tables_for_target(target_dir)
    return edge_json_path


def merge_round_results(job, workdir, round_id, reference_species, merged_target_paths):
    merged_targets = [read_json(p) for p in merged_target_paths]

    target_summary = {}
    for target in merged_targets:
        n_ok = sum(1 for r in target["results"] if r["status"] == "ok")
        n_error = sum(1 for r in target["results"] if r["status"] == "error")
        n_missing = sum(
            1 for r in target["results"]
            if r["status"] == "ok" and r.get("missing_annotations")
        )
        target_summary[target["target_species"]] = {
            "n_results": target["n_results"],
            "n_ok": n_ok,
            "n_error": n_error,
            "n_missing": n_missing,
        }

    out = {
        "round_id": round_id,
        "reference_species": reference_species,
        "targets": merged_targets,
        "target_summary": target_summary,
    }

    out_path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / "round_merged.json"
    )
    write_json(out_path, out)
    return str(out_path)


def annotate_missing_loci_and_choose_next(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    round_id,
    current_reference,
    round_merged_path,
    used_reference_species,
):
    species_list = get_species_list(species_csv)
    transcripts_by_species = load_all_transcripts(annotation_dir, annotation_suffix, species_list)
    species_loci = build_all_species_loci(transcripts_by_species)
    hal = HALAdapter(str(Path(hal_path).resolve()))

    round_merged = read_json(round_merged_path)
    missing_by_target = extract_missing_locus_payloads(round_merged)
    explained_locus_ids_by_species = collect_explained_locus_ids_from_round(round_merged)

    new_consensus_by_species = {}

    for target_species, payloads in missing_by_target.items():
        projected_transcripts = []

        for payload in payloads:
            source_species = payload["source_species"]
            source_tx_id = payload["source_transcript"]

            if source_tx_id not in transcripts_by_species[source_species]:
                continue

            seed = transcripts_by_species[source_species][source_tx_id]
            projected_exon_blocks = []

            for exon_start, exon_end in seed.exons:
                intervals = hal.project_interval(
                    source_species=seed.species,
                    target_species=target_species,
                    seqid=seed.seqid,
                    start=exon_start,
                    end=exon_end,
                    strand=seed.strand,
                    source_transcript=seed.transcript_id,
                )
                projected_exon_blocks.append(intervals)

            pts = reconstruct_projected_transcripts(seed, projected_exon_blocks)
            projected_transcripts.extend(pts)

        if not projected_transcripts:
            continue

        clusters = cluster_projected_transcripts(projected_transcripts, max_gap=0)
        consensuses = []

        for cluster in clusters:
            best_strand, _ = choose_missing_locus_strand(cluster)
            consensus = build_consensus_missing_transcript(cluster, best_strand)
            consensuses.append(consensus)

        if consensuses:
            new_consensus_by_species[target_species] = consensuses

    updated_used = append_unique_preserve_order(used_reference_species, current_reference)

    orphan_loci_by_species = {}
    pending_frontiers_by_species = {}

    for target_species in species_list:
        if target_species in updated_used:
            orphan_loci_by_species[target_species] = []
            continue

        native_loci = species_loci[target_species]
        projected_spans = collect_projected_transcript_spans_for_species(
            transcripts_by_species=transcripts_by_species,
            source_species_list=updated_used,
            target_species=target_species,
            hal=hal,
        )

        explained_ids = explained_locus_ids_by_species.get(target_species, set())
        orphan_loci = []

        for locus in native_loci:
            if locus.locus_id in explained_ids:
                continue
            if locus_overlaps_any_interval(locus, projected_spans):
                continue

            orphan_loci.append(
                {
                    "locus_id": locus.locus_id,
                    "species": locus.species,
                    "seqid": locus.seqid,
                    "start": locus.start,
                    "end": locus.end,
                    "strand": locus.strand,
                    "transcripts": locus.transcripts,
                }
            )

        orphan_loci_by_species[target_species] = orphan_loci

    for sp in species_list:
        if sp in updated_used:
            continue

        pending = []

        for c in new_consensus_by_species.get(sp, []):
            pending.append(
                {
                    "kind": "urcat_consensus",
                    "species": c.species,
                    "seqid": c.seqid,
                    "start": min(e[0] for e in c.exons),
                    "end": max(e[1] for e in c.exons),
                    "strand": c.strand,
                    "support_count": c.support_count,
                    "total_chain_score": c.total_chain_score,
                    "mean_exon_recovery": c.mean_exon_recovery,
                    "source_transcripts": c.source_transcripts,
                    "exons": c.exons,
                }
            )

        for o in orphan_loci_by_species.get(sp, []):
            pending.append(
                {
                    "kind": "orphan_native_locus",
                    **o,
                }
            )

        if pending:
            pending_frontiers_by_species[sp] = pending

    seed_species = updated_used[0]
    tree_newick = get_hal_tree_newick(hal_path)
    reference_order = compute_reference_order(
        seed_species=seed_species,
        species_tree_newick=tree_newick,
        species_list=species_list,
    )

    next_reference = pick_next_reference_from_order(
        reference_order=reference_order,
        used_reference_species=updated_used,
        pending_frontiers_by_species=pending_frontiers_by_species,
    )

    out = {
        "round_id": round_id,
        "seed_species": seed_species,
        "current_reference": current_reference,
        "reference_species": current_reference,
        "reference_order": reference_order,
        "used_reference_species": updated_used,
        "new_consensus_by_species": {
            sp: [
                {
                    "species": c.species,
                    "seqid": c.seqid,
                    "strand": c.strand,
                    "support_count": c.support_count,
                    "total_chain_score": c.total_chain_score,
                    "mean_exon_recovery": c.mean_exon_recovery,
                    "source_transcripts": c.source_transcripts,
                    "exons": c.exons,
                }
                for c in cs
            ]
            for sp, cs in new_consensus_by_species.items()
        },
        "orphan_loci_by_species": orphan_loci_by_species,
        "pending_frontiers_by_species": pending_frontiers_by_species,
        "pending_seeds_by_species": pending_frontiers_by_species,
        "next_reference_species": next_reference,
        "stop": next_reference is None,
    }

    ref_out_path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{current_reference}"
        / "post_round_decision.json"
    )
    round_out_path = (
        Path(workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / "post_round_decision.json"
    )
    write_json(ref_out_path, out)
    write_json(round_out_path, out)
    finalize_round_outputs(workdir, round_id)
    return str(ref_out_path)


def schedule_target_batches(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    round_id,
    reference_species,
    target_species,
    manifest_path,
):
    manifest = read_json(manifest_path)
    target_jobs = [j for j in manifest["jobs"] if j["target_species"] == target_species]

    batch_rvs = []
    for j in target_jobs:
        batch_job = job.addChildJobFn(
            run_project_batch,
            workdir,
            annotation_dir,
            annotation_suffix,
            hal_path,
            species_csv,
            round_id,
            reference_species,
            target_species,
            j["batch_file"],
            j["batch_id"],
            memory="4G",
            disk="4G",
        )
        batch_rvs.append(batch_job.rv())

    merge_job = job.addFollowOnJobFn(
        merge_target_results,
        workdir,
        round_id,
        reference_species,
        target_species,
        batch_rvs,
        memory="2G",
        disk="2G",
    )

    edge_job = merge_job.addFollowOnJobFn(
        run_target_edge_evidence,
        workdir,
        annotation_dir,
        annotation_suffix,
        hal_path,
        species_csv,
        merge_job.rv(),
        memory="4G",
        disk="4G",
    )

    return merge_job.rv()


def schedule_next_round(
    job,
    decision_path,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    batch_size,
):
    decision = read_json(decision_path)

    if decision["stop"]:
        return decision_path

    next_ref = decision["next_reference_species"]
    next_round_id = decision["round_id"] + 1

    frontier_job = job.addChildJobFn(
        write_species_frontier_from_pending,
        workdir,
        decision["pending_frontiers_by_species"],
        next_ref,
        next_round_id,
        memory="2G",
        disk="2G",
    )

    manifest_job = frontier_job.addFollowOnJobFn(
        write_manifest,
        workdir,
        frontier_job.rv(),
        next_ref,
        [sp for sp in get_species_list(species_csv) if sp != next_ref],
        batch_size,
        memory="2G",
        disk="2G",
    )

    next_round_job = manifest_job.addFollowOnJobFn(
        schedule_round_from_manifest,
        workdir,
        annotation_dir,
        annotation_suffix,
        hal_path,
        species_csv,
        next_ref,
        manifest_job.rv(),
        decision["used_reference_species"],
        batch_size,
        memory="2G",
        disk="2G",
    )

    return next_round_job.rv()


def schedule_round_from_manifest(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    reference_species,
    manifest_path,
    used_reference_species,
    batch_size,
):
    manifest = read_json(manifest_path)
    round_id = manifest["round_id"]

    species_list = get_species_list(species_csv)
    targets = [sp for sp in species_list if sp != reference_species]

    target_merge_rvs = []
    for target in targets:
        target_job = job.addChildJobFn(
            schedule_target_batches,
            workdir,
            annotation_dir,
            annotation_suffix,
            hal_path,
            species_csv,
            round_id,
            reference_species,
            target,
            manifest_path,
            memory="2G",
            disk="2G",
        )
        target_merge_rvs.append(target_job.rv())

    round_merge_job = job.addFollowOnJobFn(
        merge_round_results,
        workdir,
        round_id,
        reference_species,
        target_merge_rvs,
        memory="2G",
        disk="2G",
    )

    decision_job = round_merge_job.addFollowOnJobFn(
        annotate_missing_loci_and_choose_next,
        workdir,
        annotation_dir,
        annotation_suffix,
        hal_path,
        species_csv,
        round_id,
        reference_species,
        round_merge_job.rv(),
        used_reference_species,
        memory="4G",
        disk="4G",
    )

    summary_job = decision_job.addFollowOnJobFn(
        write_round_summary,
        workdir,
        round_merge_job.rv(),
        decision_job.rv(),
        memory="1G",
        disk="1G",
    )

    next_round_job = summary_job.addFollowOnJobFn(
        schedule_next_round,
        decision_job.rv(),
        workdir,
        annotation_dir,
        annotation_suffix,
        hal_path,
        species_csv,
        batch_size,
        memory="2G",
        disk="2G",
    )

    return next_round_job.rv()


def run_round_zero(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    seed_species,
    batch_size,
):
    species_list = get_species_list(species_csv)
    targets = [sp for sp in species_list if sp != seed_species]

    frontier_job = job.addChildJobFn(
        write_seed_frontier,
        workdir,
        annotation_dir,
        annotation_suffix,
        seed_species,
        memory="2G",
        disk="2G",
    )

    manifest_job = frontier_job.addFollowOnJobFn(
        write_manifest,
        workdir,
        frontier_job.rv(),
        seed_species,
        targets,
        batch_size,
        memory="2G",
        disk="2G",
    )

    round_job = manifest_job.addFollowOnJobFn(
        schedule_round_from_manifest,
        workdir,
        annotation_dir,
        annotation_suffix,
        hal_path,
        species_csv,
        seed_species,
        manifest_job.rv(),
        [seed_species],
        batch_size,
        memory="2G",
        disk="2G",
    )

    return round_job.rv()


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("--outputDir", required=True)
    parser.add_argument("--seedSpecies", required=True)
    parser.add_argument("--speciesCsv", required=True)
    parser.add_argument("--halPath", required=True)
    parser.add_argument("--annotationDir", required=True)
    parser.add_argument("--annotationSuffix", default=".test.gff3")
    parser.add_argument("--batchSize", type=int, default=200)

    options = parser.parse_args()

    output_dir = str(Path(options.outputDir).resolve())
    annotation_dir = str(Path(options.annotationDir).resolve())
    hal_path = str(Path(options.halPath).resolve())

    root = Job.wrapJobFn(
        run_round_zero,
        output_dir,
        annotation_dir,
        options.annotationSuffix,
        hal_path,
        options.speciesCsv,
        options.seedSpecies,
        options.batchSize,
        memory="2G",
        disk="2G",
    )

    with Toil(options) as toil:
        toil.start(root)

    from comparative_annotator.workflow.final_gff3 import write_all_final_species_gff3

    species_list = [x.strip() for x in options.speciesCsv.split(",") if x.strip()]

    with Toil(options) as toil:
        toil.start(root)

    write_all_final_species_gff3(
        output_dir=output_dir,
        annotation_dir=annotation_dir,
        annotation_suffix=options.annotationSuffix,
        species_list=species_list,
    )


if __name__ == "__main__":
    main()
