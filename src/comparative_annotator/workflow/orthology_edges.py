from __future__ import annotations

from pathlib import Path

from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.orthology.edge_model import (
    AnchorMap,
    Interval,
    build_and_classify_edges,
)

from comparative_annotator.workflow.sequence_prep import load_all_species_sequences

def read_json(path):
    import json

    path = Path(path)
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    import json

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def load_transcripts_for_species(annotation_dir: str, annotation_suffix: str, species: str):
    gff_path = (Path(annotation_dir) / f"{species}{annotation_suffix}").resolve()
    return load_gff3(str(gff_path), species=species)


def load_all_transcripts(annotation_dir: str, annotation_suffix: str, species_list: list[str]):
    return {
        sp: load_transcripts_for_species(annotation_dir, annotation_suffix, sp)
        for sp in species_list
    }

def canonicalize_transcript_id(tx_id: str | None) -> str | None:
    if tx_id is None:
        return None
    return tx_id.removeprefix("transcript:")


def infer_locus_id_from_transcript_id(tx_id: str | None) -> str | None:
    """
    Examples:
      transcript:Hmel202001oG10.1 -> Hmel202001oG10
      Hmel202001oG10.1            -> Hmel202001oG10
    """
    tx_id = canonicalize_transcript_id(tx_id)
    if tx_id is None:
        return None
    if "." in tx_id:
        return tx_id.rsplit(".", 1)[0]
    return tx_id

def build_all_species_loci(
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
) -> dict[str, list[SpeciesLocus]]:
    species_loci = {}
    for sp, txdict in transcripts_by_species.items():
        species_loci[sp] = build_species_loci(list(txdict.values()), species=sp)
    return species_loci


def build_transcript_lookup(
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
) -> dict[str, dict[str, CandidateTranscript]]:
    return transcripts_by_species


def build_empty_anchor_map() -> AnchorMap:
    return AnchorMap(locus_to_anchor={})


def interval_from_missing_annotation_payload(payload: dict) -> Interval | None:
    seqid = payload.get("seqid")
    start = payload.get("start")
    end = payload.get("end")
    strand = payload.get("strand")
    if seqid is None or start is None or end is None:
        return None
    return Interval(seqid=seqid, start=int(start), end=int(end), strand=strand)


def collect_candidate_pairs_from_merged_target(
    merged_target: dict,
    species_loci: dict[str, list[SpeciesLocus]],
) -> list[tuple[str, str, str, str, str, Interval | None]]:
    """
    Convert one target merged.json block into candidate_pairs expected by
    build_and_classify_edges().

    Current conservative policy:
    - use primary target locus as the main accepted candidate
    - add alternative target loci as additional candidates
    - use missing_annotations payload as projected interval hint when available

    candidate_pairs rows:
      (
        source_species,
        source_locus_id,
        target_species,
        target_locus_id,
        edge_origin,
        projected_interval,
      )
    """
    target_species = merged_target["target_species"]

    locus_ids_by_species = {
        sp: {loc.locus_id for loc in loci}
        for sp, loci in species_loci.items()
    }

    candidate_pairs = []

    for r in merged_target["results"]:
        if r.get("status") != "ok":
            continue

        source_species = r.get("source_species")
        if source_species is None:
            continue

        projected_interval = None
        
        primary = r.get("primary") or {}
        primary_target_locus = primary.get(target_species)

        if primary_target_locus and primary_target_locus in locus_ids_by_species.get(target_species, set()):
            candidate_pairs.append(
                (
                    source_species,
                    r.get("source_transcript"),
                    target_species,
                    primary_target_locus,
                    "primary",
                    projected_interval,
                )
            )

        for alt_locus_id in (r.get("alternatives") or {}).get(target_species, []):
            if alt_locus_id not in locus_ids_by_species.get(target_species, set()):
                continue
            candidate_pairs.append(
                (
                    source_species,
                    r.get("source_transcript"),
                    target_species,
                    alt_locus_id,
                    "alternative",
                    projected_interval,
                )
            )

    return candidate_pairs


def remap_source_locus_ids_from_source_transcripts(
    candidate_pairs: list[tuple[str, str, str, str, str, Interval | None]],
    merged_target: dict,
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
    species_loci: dict[str, list[SpeciesLocus]],
) -> list[tuple[str, str, str, str, str, Interval | None]]:
    """
    candidate_pairs stores source transcript IDs in the 2nd field.
    Convert them to source locus IDs using, in order:

    1. exact transcript membership in SpeciesLocus.transcripts
    2. canonicalized transcript ID (without 'transcript:' prefix)
    3. inferred locus ID from transcript ID by stripping isoform suffix
    """
    tx_to_locus: dict[tuple[str, str], str] = {}
    locus_ids_by_species: dict[str, set[str]] = {}

    for species, loci in species_loci.items():
        locus_ids_by_species[species] = {locus.locus_id for locus in loci}

        for locus in loci:
            for tx_id in locus.transcripts:
                tx_to_locus[(species, tx_id)] = locus.locus_id

                norm_tx = canonicalize_transcript_id(tx_id)
                if norm_tx is not None:
                    tx_to_locus[(species, norm_tx)] = locus.locus_id

    out = []

    for source_species, source_tx_id, target_species, target_locus_id, edge_origin, projected_interval in candidate_pairs:
        source_locus_id = None

        # 1. direct lookup
        source_locus_id = tx_to_locus.get((source_species, source_tx_id))

        # 2. canonicalized transcript ID
        if source_locus_id is None:
            norm_tx = canonicalize_transcript_id(source_tx_id)
            source_locus_id = tx_to_locus.get((source_species, norm_tx))

        # 3. infer locus ID directly from transcript ID
        if source_locus_id is None:
            inferred_locus = infer_locus_id_from_transcript_id(source_tx_id)
            if inferred_locus in locus_ids_by_species.get(source_species, set()):
                source_locus_id = inferred_locus

        if source_locus_id is None:
            continue

        out.append(
            (
                source_species,
                source_locus_id,
                target_species,
                target_locus_id,
                edge_origin,
                projected_interval,
            )
        )

    return out


def write_edge_rows_tsv(path: str, edge_rows: list[dict]):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if not edge_rows:
        with open(path, "w") as fh:
            fh.write("")
        return

    columns = sorted({k for row in edge_rows for k in row.keys()})
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in edge_rows:
            fh.write("\t".join("" if row.get(col) is None else str(row.get(col)) for col in columns) + "\n")


def build_target_edge_evidence(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_csv: str,
    merged_target_path: str,
) -> str:
    """
    Build orthology edge evidence for one target merged.json file.
    Returns path to edge_evidence.json.
    """
    species_list = [x.strip() for x in species_csv.split(",") if x.strip()]
    transcripts_by_species = load_all_transcripts(annotation_dir, annotation_suffix, species_list)
    species_loci = build_all_species_loci(transcripts_by_species)

    merged_target = read_json(merged_target_path)

    candidate_pairs = collect_candidate_pairs_from_merged_target(
        merged_target=merged_target,
        species_loci=species_loci,
    )
    candidate_pairs = remap_source_locus_ids_from_source_transcripts(
        candidate_pairs=candidate_pairs,
        merged_target=merged_target,
        transcripts_by_species=transcripts_by_species,
        species_loci=species_loci,
    )

    anchor_map = build_empty_anchor_map()

    edges, orthogroups = build_and_classify_edges(
        loci_by_species=species_loci,
        transcripts_by_species=build_transcript_lookup(transcripts_by_species),
        candidate_pairs=candidate_pairs,
        anchor_map=anchor_map,
    )

    edge_rows = [edge.to_row() for edge in edges]

    merged_target_obj = Path(merged_target_path)
    out_dir = merged_target_obj.parent
    json_path = out_dir / "edge_evidence.json"
    tsv_path = out_dir / "edge_evidence.tsv"
    orthogroups_path = out_dir / "orthogroups.json"

    payload = {
        "round_id": merged_target["round_id"],
        "reference_species": merged_target["reference_species"],
        "target_species": merged_target["target_species"],
        "n_edges": len(edge_rows),
        "edges": edge_rows,
    }

    write_json(json_path, payload)
    write_edge_rows_tsv(tsv_path, edge_rows)
    write_json(
        orthogroups_path,
        {
            "round_id": merged_target["round_id"],
            "reference_species": merged_target["reference_species"],
            "target_species": merged_target["target_species"],
            "orthogroups": {
                mode: [sorted(list(comp)) for comp in comps]
                for mode, comps in orthogroups.items()
            },
        },
    )
    return str(json_path)
