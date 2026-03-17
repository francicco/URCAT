from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import luigi


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def read_json(path: Path) -> Any:
    with open(path) as fh:
        return json.load(fh)


def write_json(path: Path, obj: Any) -> None:
    ensure_parent(path)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


class WorkflowConfig(luigi.Config):
    workdir = luigi.Parameter(default="urcat_work")
    seed_species = luigi.Parameter()
    species_csv = luigi.Parameter()
    batch_size = luigi.IntParameter(default=200)


class ParseAnnotations(luigi.Task):
    species = luigi.Parameter()

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(
            str(Path(cfg.workdir) / "parsed" / f"{self.species}.json")
        )

    def run(self):
        data = {
            "species": self.species,
            "status": "parsed",
        }
        out = Path(self.output().path)
        write_json(out, data)


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
            "from": parsed["status"],
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
        data = {
            "round_id": 0,
            "reference_species": cfg.seed_species,
            "frontier_kind": "seed_all_transcripts",
            "transcript_ids": [],
        }
        write_json(Path(self.output().path), data)


class MakeReferenceFrontier(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()

    def requires(self):
        if self.round_id == 0:
            return SeedFrontier()
        cfg = WorkflowConfig()
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

        manifest = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "batch_size": cfg.batch_size,
            "jobs": [
                {
                    "target_species": target,
                    "batch_id": 0,
                    "batch_file": f"batch_{target}_000.json",
                }
                for target in targets
            ],
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
        result = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "target_species": self.target_species,
            "batch_id": self.batch_id,
            "status": "done",
            "results": [],
        }
        write_json(Path(self.output().path), result)


class MergeTargetResults(luigi.Task):
    round_id = luigi.IntParameter()
    reference_species = luigi.Parameter()
    target_species = luigi.Parameter()

    def requires(self):
        return [
            ProjectBatch(
                round_id=self.round_id,
                reference_species=self.reference_species,
                target_species=self.target_species,
                batch_id=0,
            )
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
        merged = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "target_species": self.target_species,
            "batches": batch_results,
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
        out = {
            "round_id": self.round_id,
            "reference_species": self.reference_species,
            "targets": merged_targets,
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
