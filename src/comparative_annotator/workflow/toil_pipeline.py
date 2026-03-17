from __future__ import annotations

import json
from pathlib import Path

from toil.common import Toil
from toil.job import Job

from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus


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


def write_seed_frontier(job, workdir, annotation_dir, annotation_suffix, seed_species):
    tx = load_transcripts_for_species(annotation_dir, annotation_suffix, seed_species)
    out = {
        "round_id": 0,
        "reference_species": seed_species,
        "frontier_kind": "seed_all_transcripts",
        "transcript_ids": sorted(tx.keys()),
    }
    path = Path(workdir) / "rounds" / "round_000" / f"ref_{seed_species}" / "frontier.json"
    write_json(path, out)
    return str(path)


def write_manifest(job, workdir, frontier_path, reference_species, target_species_list, batch_size):
    frontier = read_json(frontier_path)
    transcript_ids = frontier["transcript_ids"]
    batches = list(chunked(transcript_ids, batch_size))

    jobs = []
    for target in target_species_list:
        for batch_id, batch_ids in enumerate(batches):
            batch_file = (
                Path(workdir)
                / "rounds"
                / "round_000"
                / f"ref_{reference_species}"
                / f"target_{target}"
                / f"batch_{batch_id:03d}.ids.json"
            )
            write_json(batch_file, {"transcript_ids": batch_ids})
            jobs.append(
                {
                    "target_species": target,
                    "batch_id": batch_id,
                    "batch_file": str(batch_file),
                    "n_transcripts": len(batch_ids),
                }
            )

    manifest = {
        "round_id": 0,
        "reference_species": reference_species,
        "batch_size": batch_size,
        "n_frontier_transcripts": len(transcript_ids),
        "jobs": jobs,
    }

    manifest_path = (
        Path(workdir)
        / "rounds"
        / "round_000"
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
    reference_species,
    target_species,
    batch_ids_path,
    batch_id,
):
    species_list = get_species_list(species_csv)

    transcripts_by_species = load_all_transcripts(annotation_dir, annotation_suffix, species_list)
    species_loci = build_all_species_loci(transcripts_by_species)
    hal = HALAdapter(str(Path(hal_path).resolve()))

    batch_info = read_json(batch_ids_path)
    transcript_ids = batch_info["transcript_ids"]
    ref_transcripts = transcripts_by_species[reference_species]

    results = []
    for tx_id in transcript_ids:
        if tx_id not in ref_transcripts:
            results.append(
                {
                    "source_species": reference_species,
                    "source_transcript": tx_id,
                    "target_species": target_species,
                    "status": "source_transcript_not_found",
                }
            )
            continue

        seed_transcript = ref_transcripts[tx_id]

        try:
            clocus = infer_comparative_locus(
                seed_transcript=seed_transcript,
                target_species=target_species,
                hal_adapter=hal,
                species_loci=species_loci,
                transcripts_by_species=transcripts_by_species,
            )

            result = {
                "source_species": reference_species,
                "source_transcript": tx_id,
                "target_species": target_species,
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
                "source_transcript": tx_id,
                "target_species": target_species,
                "status": "error",
                "error": str(e),
            }

        results.append(result)

    out = {
        "round_id": 0,
        "reference_species": reference_species,
        "target_species": target_species,
        "batch_id": batch_id,
        "n_transcripts": len(transcript_ids),
        "results": results,
    }

    out_path = (
        Path(workdir)
        / "rounds"
        / "round_000"
        / f"ref_{reference_species}"
        / f"target_{target_species}"
        / f"batch_{batch_id:03d}.json"
    )
    write_json(out_path, out)
    return str(out_path)


def merge_target_results(job, workdir, reference_species, target_species, batch_result_paths):
    batch_results = [read_json(p) for p in batch_result_paths]

    merged_results = []
    for batch in batch_results:
        merged_results.extend(batch["results"])

    merged = {
        "round_id": 0,
        "reference_species": reference_species,
        "target_species": target_species,
        "n_batches": len(batch_results),
        "n_results": len(merged_results),
        "results": merged_results,
    }

    out_path = (
        Path(workdir)
        / "rounds"
        / "round_000"
        / f"ref_{reference_species}"
        / f"target_{target_species}"
        / "merged.json"
    )
    write_json(out_path, merged)
    return str(out_path)


def merge_round_results(job, workdir, reference_species, merged_target_paths):
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
        "round_id": 0,
        "reference_species": reference_species,
        "targets": merged_targets,
        "target_summary": target_summary,
        "new_informative_unresolved_loci": [],
        "stop": False,
    }

    out_path = (
        Path(workdir)
        / "rounds"
        / "round_000"
        / f"ref_{reference_species}"
        / "round_merged.json"
    )
    write_json(out_path, out)
    return str(out_path)


def schedule_target_batches(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
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
        reference_species,
        target_species,
        batch_rvs,
        memory="2G",
        disk="2G",
    )
    return merge_job.rv()


def schedule_round_from_manifest(
    job,
    workdir,
    annotation_dir,
    annotation_suffix,
    hal_path,
    species_csv,
    reference_species,
    manifest_path,
):
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
            reference_species,
            target,
            manifest_path,
            memory="2G",
            disk="2G",
        )
        target_merge_rvs.append(target_job.rv())

    round_merge = job.addFollowOnJobFn(
        merge_round_results,
        workdir,
        reference_species,
        target_merge_rvs,
        memory="2G",
        disk="2G",
    )
    return round_merge.rv()


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


if __name__ == "__main__":
    main()
