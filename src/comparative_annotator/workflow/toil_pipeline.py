from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from toil.common import Toil
from toil.job import Job

from comparative_annotator.workflow.config import URCATConfig, load_urcat_config
from comparative_annotator.workflow.annotation_sources import (
    load_all_transcripts,
    build_all_species_loci,
)
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.models.transcript import CandidateTranscript
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
from comparative_annotator.workflow.workflow_paths import get_round_dir
from comparative_annotator.workflow.analysis_tables import (
    write_all_analysis_tables_for_round,
    write_all_analysis_tables_for_target,
)
from comparative_annotator.workflow.final_gff3 import write_final_species_gff3s


@dataclass
class FrontierSeedTranscript:
    transcript_id: str
    species: str
    seqid: str
    strand: str
    exons: list


# ---------------------------------------------------------------------------
# JSON / utility helpers
# ---------------------------------------------------------------------------

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


def append_unique_preserve_order(items, value):
    out = list(items)
    if value not in out:
        out.append(value)
    return out


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


# ---------------------------------------------------------------------------
# Helper functions that were missing from the user's version
# ---------------------------------------------------------------------------

def collect_explained_locus_ids_from_round(round_merged_json: dict) -> dict:
    """Return {target_species: set_of_locus_ids} already explained this round."""
    explained = {}
    for target_block in round_merged_json["targets"]:
        target_species = target_block["target_species"]
        explained_ids: set = set()
        for r in target_block["results"]:
            if r.get("status") != "ok":
                continue
            for sp, locus_id in (r.get("primary") or {}).items():
                if sp == target_species and locus_id:
                    explained_ids.add(locus_id)
            for sp, locus_ids in (r.get("alternatives") or {}).items():
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
    """Project all source transcripts into target_species; return span dicts."""
    spans = []
    for source_species in source_species_list:
        if source_species == target_species:
            continue
        for _, seed in transcripts_by_species.get(source_species, {}).items():
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
                spans.append({
                    "seqid": pt.seqid,
                    "start": pt.start,
                    "end": pt.end,
                    "strand": pt.strand,
                    "source_species": source_species,
                    "source_transcript": pt.source_transcript,
                })
    return spans


def finalize_round_outputs(output_dir: str, round_id: int) -> None:
    round_dir = get_round_dir(output_dir, round_id)
    write_all_analysis_tables_for_round(round_dir)


def write_round_summary(job, workdir, round_merged_path, decision_path):
    round_merged = read_json(round_merged_path)
    decision = read_json(decision_path)

    round_id = decision["round_id"]
    reference_species = decision["reference_species"]

    ref_dir = Path(workdir) / "rounds" / f"round_{round_id:03d}" / f"ref_{reference_species}"
    round_dir = Path(workdir) / "rounds" / f"round_{round_id:03d}"

    summary = {
        "round_id": round_id,
        "seed_species": decision["seed_species"],
        "reference_species": reference_species,
        "reference_order": decision["reference_order"],
        "used_reference_species": decision["used_reference_species"],
        "next_reference_species": decision["next_reference_species"],
        "stop": decision["stop"],
    }
    for json_path in [ref_dir / "summary.json", round_dir / "summary.json"]:
        write_json(json_path, summary)

    header = [
        "round_id", "seed_species", "reference_species", "target_species",
        "n_processed", "n_ok", "n_error", "n_primary", "n_alternative_only",
        "n_missing", "n_strand_conflict", "n_new_consensus", "n_orphan_loci",
        "n_pending_frontier", "next_reference_species", "stop",
    ]

    targets_summary = {}
    for target_block in round_merged.get("targets", []):
        sp = target_block["target_species"]
        results = target_block["results"]
        targets_summary[sp] = {
            "n_processed": len(results),
            "n_ok": sum(1 for r in results if r.get("status") == "ok"),
            "n_error": sum(1 for r in results if r.get("status") == "error"),
            "n_primary": sum(1 for r in results if r.get("status") == "ok" and r.get("primary")),
            "n_alternative_only": sum(
                1 for r in results
                if r.get("status") == "ok" and not r.get("primary") and r.get("alternatives")
            ),
            "n_missing": sum(
                1 for r in results if r.get("status") == "ok" and r.get("missing_annotations")
            ),
            "n_strand_conflict": sum(
                1 for r in results if r.get("status") == "ok" and r.get("strand_conflicts")
            ),
            "n_new_consensus": len(decision.get("new_consensus_by_species", {}).get(sp, [])),
            "n_orphan_loci": len(decision.get("orphan_loci_by_species", {}).get(sp, [])),
            "n_pending_frontier": len(decision.get("pending_frontiers_by_species", {}).get(sp, [])),
        }

    for tsv_path in [ref_dir / "summary.tsv", round_dir / "summary.tsv"]:
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(tsv_path, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for sp, stats in targets_summary.items():
                fh.write("\t".join([
                    str(round_id), str(decision["seed_species"]),
                    str(reference_species), str(sp),
                    str(stats["n_processed"]), str(stats["n_ok"]),
                    str(stats["n_error"]), str(stats["n_primary"]),
                    str(stats["n_alternative_only"]), str(stats["n_missing"]),
                    str(stats["n_strand_conflict"]), str(stats["n_new_consensus"]),
                    str(stats["n_orphan_loci"]), str(stats["n_pending_frontier"]),
                    str(decision["next_reference_species"]), str(decision["stop"]),
                ]) + "\n")

    return str(ref_dir / "summary.json")


# ---------------------------------------------------------------------------
# Adapter for build_target_edge_evidence
# orthology_edges.build_target_edge_evidence still uses annotation_dir/suffix/
# species_csv. This wrapper converts cfg.annotation_paths to those args.
# ---------------------------------------------------------------------------

def _run_target_edge_evidence_impl(workdir: str, cfg: URCATConfig, merged_target_path: str) -> str:
    import comparative_annotator.workflow.orthology_edges as oe
    _orig = oe.load_all_transcripts

    def _patched(annotation_dir, annotation_suffix, species_list):
        return load_all_transcripts(cfg.annotation_paths, species_list)

    oe.load_all_transcripts = _patched
    try:
        edge_json_path = oe.build_target_edge_evidence(
            workdir=workdir,
            annotation_dir="",
            annotation_suffix="",
            hal_path=cfg.hal_path,
            species_csv=",".join(cfg.species_list),
            merged_target_path=merged_target_path,
        )
    finally:
        oe.load_all_transcripts = _orig

    return edge_json_path


# ---------------------------------------------------------------------------
# Toil job functions
# ---------------------------------------------------------------------------

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


def write_seed_frontier(job, workdir, annotation_paths, seed_species):
    gff_path = annotation_paths.get(seed_species, "")
    tx = load_gff3(gff_path, species=seed_species) if gff_path else {}
    out = {
        "round_id": 0,
        "reference_species": seed_species,
        "frontier_kind": "seed_all_transcripts",
        "transcript_ids": sorted(tx.keys()),
    }
    path = Path(workdir) / "rounds" / "round_000" / f"ref_{seed_species}" / "frontier.json"
    write_json(path, out)
    return str(path)


def write_species_frontier_from_pending(
    job, workdir, pending_frontiers_by_species, reference_species, round_id,
):
    pending = pending_frontiers_by_species.get(reference_species, [])
    frontier = {
        "round_id": round_id,
        "reference_species": reference_species,
        "frontier_kind": "mixed_pending",
        "items": pending,
    }
    path = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}"
        / f"ref_{reference_species}" / "frontier.json"
    )
    write_json(path, frontier)
    return str(path)


def write_manifest(job, workdir, frontier_path, reference_species, target_species_list, batch_size):
    frontier = read_json(frontier_path)
    round_id = frontier["round_id"]

    if frontier["frontier_kind"] == "seed_all_transcripts":
        items = [{"kind": "native_transcript", "transcript_id": tid}
                 for tid in frontier["transcript_ids"]]
    elif frontier["frontier_kind"] == "mixed_pending":
        items = frontier["items"]
    else:
        raise ValueError(f"Unknown frontier kind: {frontier['frontier_kind']}")

    batches = list(chunked(items, batch_size))
    jobs = []
    for target in target_species_list:
        for batch_id, batch_items in enumerate(batches):
            batch_file = (
                Path(workdir) / "rounds" / f"round_{round_id:03d}"
                / f"ref_{reference_species}" / f"target_{target}"
                / f"batch_{batch_id:03d}.items.json"
            )
            write_json(batch_file, {"items": batch_items})
            jobs.append({
                "round_id": round_id,
                "target_species": target,
                "batch_id": batch_id,
                "batch_file": str(batch_file),
                "n_items": len(batch_items),
            })

    manifest = {
        "round_id": round_id,
        "reference_species": reference_species,
        "batch_size": batch_size,
        "jobs": jobs,
    }
    manifest_path = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}"
        / f"ref_{reference_species}" / "manifest.json"
    )
    write_json(manifest_path, manifest)
    return str(manifest_path)


def run_project_batch(
    job, workdir, cfg: URCATConfig, round_id,
    reference_species, target_species, batch_items_path, batch_id,
):
    transcripts_by_species = load_all_transcripts(cfg.annotation_paths, cfg.species_list)
    species_loci = build_all_species_loci(transcripts_by_species)
    hal = HALAdapter(str(Path(cfg.hal_path).resolve()))

    items = read_json(batch_items_path)["items"]
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
        Path(workdir) / "rounds" / f"round_{round_id:03d}"
        / f"ref_{reference_species}" / f"target_{target_species}"
        / f"batch_{batch_id:03d}.json"
    )
    write_json(out_path, out)
    return str(out_path)


def merge_target_results(
    job, workdir, round_id, reference_species, target_species, batch_result_paths,
):
    merged_results = []
    for p in batch_result_paths:
        merged_results.extend(read_json(p)["results"])

    merged = {
        "round_id": round_id,
        "reference_species": reference_species,
        "target_species": target_species,
        "n_batches": len(batch_result_paths),
        "n_results": len(merged_results),
        "results": merged_results,
    }
    out_path = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}"
        / f"ref_{reference_species}" / f"target_{target_species}" / "merged.json"
    )
    write_json(out_path, merged)
    return str(out_path)


def run_target_edge_evidence(job, workdir: str, cfg: URCATConfig, merged_target_path: str) -> str:
    edge_json_path = _run_target_edge_evidence_impl(workdir, cfg, merged_target_path)
    write_all_analysis_tables_for_target(Path(edge_json_path).parent)
    return edge_json_path


def merge_round_results(job, workdir, round_id, reference_species, merged_target_paths):
    merged_targets = [read_json(p) for p in merged_target_paths]
    target_summary = {}
    for target in merged_targets:
        results = target["results"]
        target_summary[target["target_species"]] = {
            "n_results": target["n_results"],
            "n_ok": sum(1 for r in results if r["status"] == "ok"),
            "n_error": sum(1 for r in results if r["status"] == "error"),
            "n_missing": sum(
                1 for r in results if r["status"] == "ok" and r.get("missing_annotations")
            ),
        }

    out = {
        "round_id": round_id,
        "reference_species": reference_species,
        "targets": merged_targets,
        "target_summary": target_summary,
    }
    out_path = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}"
        / f"ref_{reference_species}" / "round_merged.json"
    )
    write_json(out_path, out)
    return str(out_path)


def annotate_missing_loci_and_choose_next(
    job, workdir, cfg: URCATConfig, round_id,
    current_reference, round_merged_path, used_reference_species,
):
    from comparative_annotator.workflow.new_loci_gff3 import write_new_loci_gff3
    from comparative_annotator.workflow.fragmented_loci_table import write_fragmented_loci_table
    from comparative_annotator.workflow.fragmented_projection import summarize_projected_blocks

    species_list = cfg.species_list
    transcripts_by_species = load_all_transcripts(cfg.annotation_paths, species_list)
    species_loci = build_all_species_loci(transcripts_by_species)
    hal = HALAdapter(str(Path(cfg.hal_path).resolve()))

    round_merged = read_json(round_merged_path)
    missing_by_target = extract_missing_locus_payloads(round_merged)
    explained_locus_ids_by_species = collect_explained_locus_ids_from_round(round_merged)

    new_consensus_by_species = {}
    fragmented_by_species = {}
    skipped_missing_payloads = []

    for target_species, payloads in missing_by_target.items():
        projected_transcripts = []
        projected_blocks_for_summary = []

        for payload in payloads:
            source_species = payload["source_species"]
            source_tx_id = payload["source_transcript"]

            if source_tx_id not in transcripts_by_species.get(source_species, {}):
                skipped_missing_payloads.append({
                    "target_species": target_species,
                    "source_species": source_species,
                    "source_transcript": source_tx_id,
                    "reason": "source_transcript_not_found",
                })
                continue

            seed = transcripts_by_species[source_species][source_tx_id]
            projected_exon_blocks = []

            for exon_idx, (exon_start, exon_end) in enumerate(seed.exons, start=1):
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

                for iv in intervals:
                    projected_blocks_for_summary.append({
                        "source_species": seed.species,
                        "source_transcript": seed.transcript_id,
                        "target_species": target_species,
                        "source_seqid": seed.seqid,
                        "source_strand": seed.strand,
                        "source_exon_number": exon_idx,
                        "target_seqid": iv.seqid,
                        "target_start": iv.start,
                        "target_end": iv.end,
                        "target_strand": iv.strand,
                        "chain_score": getattr(iv, "chain_score", None),
                    })

            pts = reconstruct_projected_transcripts(seed, projected_exon_blocks)
            projected_transcripts.extend(pts)

        if projected_transcripts:
            clusters = cluster_projected_transcripts(projected_transcripts, max_gap=0)
            consensuses = [
                build_consensus_missing_transcript(cluster, choose_missing_locus_strand(cluster)[0])
                for cluster in clusters
            ]
            consensuses = [c for c in consensuses if c is not None]
            if consensuses:
                new_consensus_by_species[target_species] = consensuses

        grouped_fragmented = {}
        for row in projected_blocks_for_summary:
            key = (row["source_species"], tuple(sorted([row["source_transcript"]])))
            grouped_fragmented.setdefault(key, []).append(row)

        fragmented_by_species[target_species] = [
            summarize_projected_blocks(
                source_species=src_sp,
                target_species=target_species,
                source_transcripts=list(src_txs),
                blocks=blocks,
            )
            for (src_sp, src_txs), blocks in grouped_fragmented.items()
        ]

    updated_used = append_unique_preserve_order(used_reference_species, current_reference)
    orphan_loci_by_species = {}
    pending_frontiers_by_species = {}

    for target_species in species_list:
        if target_species in updated_used:
            orphan_loci_by_species[target_species] = []
            continue

        native_loci = species_loci.get(target_species, [])
        projected_spans = collect_projected_transcript_spans_for_species(
            transcripts_by_species=transcripts_by_species,
            source_species_list=updated_used,
            target_species=target_species,
            hal=hal,
        )
        explained_ids = explained_locus_ids_by_species.get(target_species, set())

        orphan_loci_by_species[target_species] = [
            {
                "locus_id": locus.locus_id,
                "species": locus.species,
                "seqid": locus.seqid,
                "start": locus.start,
                "end": locus.end,
                "strand": locus.strand,
                "transcripts": locus.transcripts,
            }
            for locus in native_loci
            if locus.locus_id not in explained_ids
            and not locus_overlaps_any_interval(locus, projected_spans)
        ]

    for sp in species_list:
        if sp in updated_used:
            continue
        pending = [
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
            for c in new_consensus_by_species.get(sp, [])
        ] + [
            {"kind": "orphan_native_locus", **o}
            for o in orphan_loci_by_species.get(sp, [])
        ]
        if pending:
            pending_frontiers_by_species[sp] = pending

    seed_species = updated_used[0]
    tree_newick = get_hal_tree_newick(cfg.hal_path)
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
                    "species": c.species, "seqid": c.seqid, "strand": c.strand,
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
        "fragmented_by_species": fragmented_by_species,
        "orphan_loci_by_species": orphan_loci_by_species,
        "pending_frontiers_by_species": pending_frontiers_by_species,
        "pending_seeds_by_species": pending_frontiers_by_species,
        "skipped_missing_payloads": skipped_missing_payloads,
        "next_reference_species": next_reference,
        "stop": next_reference is None,
    }

    ref_out_path = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}"
        / f"ref_{current_reference}" / "post_round_decision.json"
    )
    round_out_path = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}" / "post_round_decision.json"
    )
    write_json(ref_out_path, out)
    write_json(round_out_path, out)

    round_ref_dir = (
        Path(workdir) / "rounds" / f"round_{round_id:03d}" / f"ref_{current_reference}"
    )
    round_ref_dir.mkdir(parents=True, exist_ok=True)

    for target_species, loci in new_consensus_by_species.items():
        from comparative_annotator.workflow.new_loci_gff3 import write_new_loci_gff3
        write_new_loci_gff3(str(round_ref_dir / f"{target_species}.new_loci.gff3"), target_species, loci)

    for target_species, loci in fragmented_by_species.items():
        write_fragmented_loci_table(
            str(round_ref_dir / f"{target_species}.fragmented_loci.tsv"),
            round_id, current_reference, target_species, loci,
        )

    finalize_round_outputs(workdir, round_id)
    return str(ref_out_path)


def schedule_target_batches(
    job, workdir, cfg: URCATConfig, round_id,
    reference_species, target_species, manifest_path,
):
    manifest = read_json(manifest_path)
    target_jobs = [j for j in manifest["jobs"] if j["target_species"] == target_species]

    batch_rvs = []
    for j in target_jobs:
        batch_job = job.addChildJobFn(
            run_project_batch,
            workdir, cfg, round_id, reference_species, target_species,
            j["batch_file"], j["batch_id"],
            memory="4G", disk="4G",
        )
        batch_rvs.append(batch_job.rv())

    merge_job = job.addFollowOnJobFn(
        merge_target_results,
        workdir, round_id, reference_species, target_species, batch_rvs,
        memory="2G", disk="2G",
    )
    merge_job.addFollowOnJobFn(
        run_target_edge_evidence,
        workdir, cfg, merge_job.rv(),
        memory="4G", disk="4G",
    )
    return merge_job.rv()


def schedule_next_round(job, decision_path, workdir, cfg: URCATConfig):
    decision = read_json(decision_path)
    if decision["stop"]:
        return decision_path

    next_ref = decision["next_reference_species"]
    next_round_id = decision["round_id"] + 1

    frontier_job = job.addChildJobFn(
        write_species_frontier_from_pending,
        workdir, decision["pending_frontiers_by_species"], next_ref, next_round_id,
        memory="2G", disk="2G",
    )
    manifest_job = frontier_job.addFollowOnJobFn(
        write_manifest,
        workdir, frontier_job.rv(), next_ref,
        [sp for sp in cfg.species_list if sp != next_ref],
        cfg.batch_size,
        memory="2G", disk="2G",
    )
    next_round_job = manifest_job.addFollowOnJobFn(
        schedule_round_from_manifest,
        workdir, cfg, next_ref, manifest_job.rv(),
        decision["used_reference_species"],
        memory="2G", disk="2G",
    )
    return next_round_job.rv()


def schedule_round_from_manifest(
    job, workdir, cfg: URCATConfig, reference_species, manifest_path, used_reference_species,
):
    manifest = read_json(manifest_path)
    round_id = manifest["round_id"]
    targets = [sp for sp in cfg.species_list if sp != reference_species]

    target_merge_rvs = []
    for target in targets:
        target_job = job.addChildJobFn(
            schedule_target_batches,
            workdir, cfg, round_id, reference_species, target, manifest_path,
            memory="2G", disk="2G",
        )
        target_merge_rvs.append(target_job.rv())

    round_merge_job = job.addFollowOnJobFn(
        merge_round_results,
        workdir, round_id, reference_species, target_merge_rvs,
        memory="2G", disk="2G",
    )
    decision_job = round_merge_job.addFollowOnJobFn(
        annotate_missing_loci_and_choose_next,
        workdir, cfg, round_id, reference_species,
        round_merge_job.rv(), used_reference_species,
        memory="4G", disk="4G",
    )
    summary_job = decision_job.addFollowOnJobFn(
        write_round_summary,
        workdir, round_merge_job.rv(), decision_job.rv(),
        memory="1G", disk="1G",
    )
    next_round_job = summary_job.addFollowOnJobFn(
        schedule_next_round,
        decision_job.rv(), workdir, cfg,
        memory="2G", disk="2G",
    )
    return next_round_job.rv()


def run_round_zero(job, workdir: str, cfg: URCATConfig):
    """Schedule round 0: project all seed-species transcripts into every target."""
    seed_species = cfg.seed_species
    targets = [sp for sp in cfg.species_list if sp != seed_species]

    frontier_job = job.addChildJobFn(
        write_seed_frontier,
        workdir, cfg.annotation_paths, seed_species,
        memory="2G", disk="2G",
    )
    manifest_job = frontier_job.addFollowOnJobFn(
        write_manifest,
        workdir, frontier_job.rv(), seed_species, targets, cfg.batch_size,
        memory="2G", disk="2G",
    )
    round_job = manifest_job.addFollowOnJobFn(
        schedule_round_from_manifest,
        workdir, cfg, seed_species, manifest_job.rv(), [seed_species],
        memory="2G", disk="2G",
    )
    return round_job.rv()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # getDefaultArgumentParser() uses configargparse and already registers
    # --config for INI-file support. Use --urcatConfig for URCAT's own config
    # to avoid the conflict we hit with --config.
    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--urcatConfig", required=True, help="Path to URCAT INI config file")
    parser.add_argument("--outputDir", required=True, help="Workflow output directory")

    args = parser.parse_args()

    cfg = load_urcat_config(args.urcatConfig)
    output_dir = str(Path(args.outputDir).resolve())

    root = Job.wrapJobFn(
        run_round_zero,
        output_dir,
        cfg,
        memory="2G",
        disk="2G",
    )

    with Toil(args) as toil:
        toil.start(root)

    write_final_species_gff3s(
        output_dir=output_dir,
        annotation_paths=cfg.annotation_paths,
        species_list=cfg.species_list,
    )


if __name__ == "__main__":
    main()
