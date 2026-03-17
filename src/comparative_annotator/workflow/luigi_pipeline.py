from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import luigi

from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def read_json(path: Path) -> Any:
    with open(path) as fh:
        return json.load(fh)


def write_json(path: Path, obj: Any) -> None:
    ensure_parent(path)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def chunked(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


class WorkflowConfig(luigi.Config):
    workdir = luigi.Parameter(default="urcat_work")
    seed_species = luigi.Parameter()
    species_csv = luigi.Parameter()
    batch_size = luigi.IntParameter(default=200)

    hal_path = luigi.Parameter()
    annotation_dir = luigi.Parameter(default="data")

    annotation_suffix = luigi.Parameter(default=".test.gff3")


def load_transcripts_for_species(species: str):
    cfg = WorkflowConfig()
    gff_path = Path(cfg.annotation_dir) / f"{species}{cfg.annotation_suffix}"
    return load_gff3(str(gff_path), species=species)


def load_all_transcripts(species_list):
    return {sp: load_transcripts_for_species(sp) for sp in species_list}


def build_all_species_loci(transcripts_by_species):
    species_loci = {}
    for sp, txdict in transcripts_by_species.items():
        species_loci[sp] = build_species_loci(list(txdict.values()), species=sp)
    return species_loci


class ParseAnnotations(luigi.Task):
    species = luigi.Parameter()

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(Path(cfg.workdir) / "parsed" / f"{self.species}.json")
        )

    def run(self):
        tx = load_transcripts_for_species(self.species)
        data = {
            "species": self.species,
            "status": "parsed",
            "n_transcripts": len(tx),
            "transcript_ids": sorted(tx.keys()),
        }
        write_json(Path(self.output().path), data)


class BuildSpeciesLoci(luigi.Task):
    species = luigi.Parameter()

    def requires(self):
        return ParseAnnotations(species=self.species)

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(Path(cfg.workdir) / "loci" / f"{self.species}.json")
        )

    def run(self):
        parsed = read_json(Path(self.input().path))
        data = {
            "species": self.species,
            "status": "loci_built",
            "n_transcripts": parsed["n_transcripts"],
        }
        write_json(Path(self.output().path), data)


class PlanTraversal(luigi.Task):
    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(str(Path(cfg.workdir) / "plan" / "traversal.json"))

    def run(self):
        cfg = WorkflowConfig()
        species = [x.strip() for x in cfg.species_csv.split(",") if x.strip()]
        plan = {
            "seed_species": cfg.seed_species,
            "species_order": species,
        }
        write_json(Path(self.output().path), plan)


class SeedFrontier(luigi.Task):
    def requires(self):
        cfg = WorkflowConfig()
        return {
            "plan": PlanTraversal(),
            "parsed": ParseAnnotations(species=cfg.seed_species),
        }

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(str(Path(cfg.workdir) / "rounds" / "round_000" / "frontier.json"))

    def run(self):
        cfg = WorkflowConfig()
        parsed = read_json(Path(self.input()["parsed"].path))
        data = {
            "round_id": 0,
            "reference_species": cfg.seed_species,
            "frontier_kind": "seed_all_transcripts",
            "transcript_ids": parsed["transcript_ids"],
        }
        write_json(Path(self.output().path), data)


class MakeReferenceFrontier(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()

    def requires(self):
        if self.round_id == 0:
            return SeedFrontier()
        return MergeRoundResults(
            round_id=self.round_id - 1,
            reference_species=self.reference_species,
        )

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(
                Path(cfg.workdir)
                / "rounds"
                / f"round_{self.round_id:03d}"
                / f"ref_{self.reference_species}"
                / "frontier.json"
            )
        )

    def run(self):
        if self.round_id == 0:
            seed = read_json(Path(self.input().path))
            data = {
                "round_id": self.round_id,
                "reference_species": self.reference_species,
                "frontier_kind": seed["frontier_kind"],
                "transcript_ids": seed["transcript_ids"],
            }
        else:
            data = {
                "round_id": self.round_id,
                "reference_species": self.reference_species,
                "frontier_kind": "unresolved_or_new",
                "transcript_ids": [],
            }

        write_json(Path(self.output().path), data)


class MakeBatchManifest(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()

    def requires(self):
        cfg = WorkflowConfig()
        species = [x.strip() for x in cfg.species_csv.split(",") if x.strip()]
        reqs = {"frontier": MakeReferenceFrontier(round_id=self.round_id, reference_species=self.reference_species)}
        reqs.update({f"loci_{sp}": BuildSpeciesLoci(species=sp) for sp in species})
        return reqs

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(
                Path(cfg.workdir)
                / "rounds"
                / f"round_{self.round_id:03d}"
                / f"ref_{self.reference_species}"
                / "manifest.json"
            )
        )

    def run(self):
        cfg = WorkflowConfig()
        species = [x.strip() for x in cfg.species_csv.split(",") if x.strip()]
        targets = [sp for sp in species if sp != self.reference_species]

        frontier = read_json(Path(self.input()["frontier"].path))
        transcript_ids = frontier["transcript_ids"]
        batches = list(chunked(transcript_ids, cfg.batch_size))

        jobs = []
        for target in targets:
            for batch_id, batch_ids in enumerate(batches):
                batch_file = (
                    Path(cfg.workdir)
                    / "rounds"
                    / f"round_{self.round_id:03d}"
                    / f"ref_{self.reference_species}"
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
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "batch_size": cfg.batch_size,
            "n_frontier_transcripts": len(transcript_ids),
            "jobs": jobs,
        }
        write_json(Path(self.output().path), manifest)


class ProjectBatch(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()
    target_species = luigi.Parameter()
    batch_id = luigi.IntParameter()

    def requires(self):
        return MakeBatchManifest(
            round_id=self.round_id,
            reference_species=self.reference_species,
        )

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(
                Path(cfg.workdir)
                / "rounds"
                / f"round_{self.round_id:03d}"
                / f"ref_{self.reference_species}"
                / f"target_{self.target_species}"
                / f"batch_{self.batch_id:03d}.json"
            )
        )

    def run(self):
        cfg = WorkflowConfig()
        species = [x.strip() for x in cfg.species_csv.split(",") if x.strip()]

        transcripts_by_species = load_all_transcripts(species)
        species_loci = build_all_species_loci(transcripts_by_species)
        hal = HALAdapter(cfg.hal_path)

        batch_ids_path = (
            Path(cfg.workdir)
            / "rounds"
            / f"round_{self.round_id:03d}"
            / f"ref_{self.reference_species}"
            / f"target_{self.target_species}"
            / f"batch_{self.batch_id:03d}.ids.json"
        )
        batch_info = read_json(batch_ids_path)
        transcript_ids = batch_info["transcript_ids"]

        ref_transcripts = transcripts_by_species[self.reference_species]

        results = []
        for tx_id in transcript_ids:
            if tx_id not in ref_transcripts:
                results.append(
                    {
                        "source_species": self.reference_species,
                        "source_transcript": tx_id,
                        "target_species": self.target_species,
                        "status": "source_transcript_not_found",
                    }
                )
                continue

            seed_transcript = ref_transcripts[tx_id]

            try:
                clocus = infer_comparative_locus(
                    seed_transcript=seed_transcript,
                    target_species=self.target_species,
                    hal_adapter=hal,
                    species_loci=species_loci,
                    transcripts_by_species=transcripts_by_species,
                )

                result = {
                    "source_species": self.reference_species,
                    "source_transcript": tx_id,
                    "target_species": self.target_species,
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
                    "source_species": self.reference_species,
                    "source_transcript": tx_id,
                    "target_species": self.target_species,
                    "status": "error",
                    "error": str(e),
                }

            results.append(result)

        out = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "target_species": self.target_species,
            "batch_id": self.batch_id,
            "n_transcripts": len(transcript_ids),
            "results": results,
        }
        write_json(Path(self.output().path), out)


class MergeTargetResults(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()
    target_species = luigi.Parameter()

    def requires(self):
        manifest = read_json(
            Path(
                MakeBatchManifest(
                    round_id=self.round_id,
                    reference_species=self.reference_species,
                ).output().path
            )
        )
        target_jobs = [
            j for j in manifest["jobs"]
            if j["target_species"] == self.target_species
        ]
        return [
            ProjectBatch(
                round_id=self.round_id,
                reference_species=self.reference_species,
                target_species=self.target_species,
                batch_id=j["batch_id"],
            )
            for j in target_jobs
        ]

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(
                Path(cfg.workdir)
                / "rounds"
                / f"round_{self.round_id:03d}"
                / f"ref_{self.reference_species}"
                / f"target_{self.target_species}"
                / "merged.json"
            )
        )

    def run(self):
        batch_results = [read_json(Path(inp.path)) for inp in self.input()]

        merged_results = []
        for batch in batch_results:
            merged_results.extend(batch["results"])

        merged = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "target_species": self.target_species,
            "n_batches": len(batch_results),
            "n_results": len(merged_results),
            "results": merged_results,
        }
        write_json(Path(self.output().path), merged)


class MergeRoundResults(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()

    def requires(self):
        cfg = WorkflowConfig()
        species = [x.strip() for x in cfg.species_csv.split(",") if x.strip()]
        targets = [sp for sp in species if sp != self.reference_species]
        return [
            MergeTargetResults(
                round_id=self.round_id,
                reference_species=self.reference_species,
                target_species=target,
            )
            for target in targets
        ]

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(
                Path(cfg.workdir)
                / "rounds"
                / f"round_{self.round_id:03d}"
                / f"ref_{self.reference_species}"
                / "round_merged.json"
            )
        )

    def run(self):
        merged_targets = [read_json(Path(inp.path)) for inp in self.input()]

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
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "targets": merged_targets,
            "target_summary": target_summary,
            "new_informative_unresolved_loci": [],
            "stop": False,
        }
        write_json(Path(self.output().path), out)


class RunProgressiveRound(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()

    def requires(self):
        return MergeRoundResults(
            round_id=self.round_id,
            reference_species=self.reference_species,
        )

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(
                Path(cfg.workdir)
                / "rounds"
                / f"round_{self.round_id:03d}"
                / f"ref_{self.reference_species}"
                / "complete.json"
            )
        )

    def run(self):
        merged = read_json(Path(self.input().path))
        data = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "status": "complete",
            "stop": merged["stop"],
            "target_summary": merged["target_summary"],
        }
        write_json(Path(self.output().path), data)


class RunInitialWorkflow(luigi.WrapperTask):
    def requires(self):
        cfg = WorkflowConfig()
        species = [x.strip() for x in cfg.species_csv.split(",") if x.strip()]
        yield PlanTraversal()
        for sp in species:
            yield BuildSpeciesLoci(species=sp)
        yield RunProgressiveRound(round_id=0, reference_species=cfg.seed_species)
