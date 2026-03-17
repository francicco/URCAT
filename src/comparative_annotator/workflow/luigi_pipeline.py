from __future__ import annotations

import json
import math
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


def get_species_list():
    cfg = WorkflowConfig()
    return [x.strip() for x in cfg.species_csv.split(",") if x.strip()]


def get_frontier_path(round_id: int, reference_species: str) -> Path:
    cfg = WorkflowConfig()
    return (
        Path(cfg.workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / "frontier.json"
    )


def get_manifest_path(round_id: int, reference_species: str) -> Path:
    cfg = WorkflowConfig()
    return (
        Path(cfg.workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / "manifest.json"
    )


def get_batch_ids_path(round_id: int, reference_species: str, target_species: str, batch_id: int) -> Path:
    cfg = WorkflowConfig()
    return (
        Path(cfg.workdir)
        / "rounds"
        / f"round_{round_id:03d}"
        / f"ref_{reference_species}"
        / f"target_{target_species}"
        / f"batch_{batch_id:03d}.ids.json"
    )


def get_number_of_batches(round_id: int, reference_species: str) -> int:
    cfg = WorkflowConfig()
    frontier_path = get_frontier_path(round_id, reference_species)
    frontier = read_json(frontier_path)
    n = len(frontier["transcript_ids"])
    if n == 0:
        return 0
    return math.ceil(n / cfg.batch_size)


class ParseAnnotations(luigi.Task):
    species = luigi.Parameter()

    def output(self):
        cfg = WorkflowConfig()
        return luigi.LocalTarget(str(Path(cfg.workdir) / "parsed" / f"{self.species}.json"))

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
        return luigi.LocalTarget(str(Path(cfg.workdir) / "loci" / f"{self.species}.json"))

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
        species = get_species_list()
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
        return luigi.LocalTarget(str(get_frontier_path(self.round_id, self.reference_species)))

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
        species = get_species_list()
        reqs = {"frontier": MakeReferenceFrontier(round_id=self.round_id, reference_species=self.reference_species)}
        reqs.update({f"loci_{sp}": BuildSpeciesLoci(species=sp) for sp in species})
        return reqs

    def output(self):
        return luigi.LocalTarget(str(get_manifest_path(self.round_id, self.reference_species)))

    def run(self):
        cfg = WorkflowConfig()
        species = get_species_list()
        targets = [sp for sp in species if sp != self.reference_species]

        frontier = read_json(Path(self.input()["frontier"].path))
        transcript_ids = frontier["transcript_ids"]
        batches = list(chunked(transcript_ids, cfg.batch_size))

        jobs = []
        for target in targets:
            for batch_id, batch_ids in enumerate(batches):
                batch_file = get_batch_ids_path(self.round_id, self.reference_species, target, batch_id)
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
       
