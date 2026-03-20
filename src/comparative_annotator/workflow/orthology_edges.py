from __future__ import annotations

import json
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
from comparative_annotator.workflow.sequence_prep import (
    load_all_species_sequences,
    prepare_diamond_inputs,
    run_diamond,
    load_diamond_results,
)


def read_json(path):
    path = Path(path)
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


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
            fh.write(
                "\t".join("" if row.get(col) is None else str(row.get(col)) for col in columns)
                + "\n"
            )


def load_transcripts_for_species(annotation_dir: str, annotation_suffix: str, species: str):
    gff_path = (Path(annotation_dir) / f"{species}{annotation_suffix}").resolve()
    return load_gff3(str(gff_path), species=species)


def load_all_transcripts(annotation_dir: str, annotation_suffix: str, species_list: list[str]):
    return {
        sp: load_transcripts_for_species(annotation_dir, annotation_suffix, sp)
        for sp in species_list
    }


def build_all_species_loci(
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
) -> dict[str, list[SpeciesLocus]]:
    return {
        sp: build_species_loci(list(txdict.values()), species=sp)
        for sp, txdict in transcripts_by_species.items()
    }


def build_transcript_lookup(
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
) -> dict[str, dict[str, CandidateTranscript]]:
    return transcripts_by_species


def build_empty_anchor_map() -> AnchorMap:
    return AnchorMap(locus_to_anchor={})


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


def parse_missing_annotation_token(token: str | None) -> Interval | None:
    """
    Parse strings like:
      Diul2100:6310-6731:+
      Eisa2100:1492633-1492980:+
    """
    if token is None:
        return None

    try:
        seqid, coords, strand = token.split(":")
        start_s, end_s = coords.split("-")
        return Interval(
            seqid=seqid,
            start=int(start_s),
            end=int(end_s),
            strand=strand,
        )
    except Exception:
        return None


def first_missing_interval_for_target(result_row: dict, target_species: str) -> Interval | None:
    missing = result_row.get("missing_annotations") or {}
    tokens = missing.get(target_species) or []
    if not tokens:
        return None
    return parse_missing_annotation_token(tokens[0])


def collect_candidate_pairs_from_merged_target(
    merged_target: dict,
    species_loci: dict[str, list[SpeciesLocus]],
) -> list[tuple[str, str, str, str, str, Interval | None]]:
    """
    Convert one target merged.json block into candidate_pairs expected by
    build_and_classify_edges().

    candidate_pairs rows:
      (
        source_species,
        source_transcript_id,
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
        source_transcript = r.get("source_transcript")
        if source_species is None or source_transcript is None:
            continue

        projected_interval = first_missing_interval_for_target(r, target_species)

        primary = r.get("primary") or {}
        primary_target_locus = primary.get(target_species)

        if (
            primary_target_locus is not None
            and primary_target_locus in locus_ids_by_species.get(target_species, set())
        ):
            candidate_pairs.append(
                (
                    source_species,
                    source_transcript,
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
                    source_transcript,
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

        source_locus_id = tx_to_locus.get((source_species, source_tx_id))

        if source_locus_id is None:
            norm_tx = canonicalize_transcript_id(source_tx_id)
            source_locus_id = tx_to_locus.get((source_species, norm_tx))

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


def build_diamond_cache_for_target(
    workdir: str,
    target_species: str,
    source_species_list: list[str],
    sequences_by_species: dict[str, dict[str, dict[str, str]]],
) -> dict[tuple[str, str], dict[tuple[str, str], dict]]:
    """
    Build/load DIAMOND hits for each source->target pair.

    Returns:
      {
        (source_species, target_species): {
          (qseqid, sseqid): {"pid": ..., "aln_len": ..., "bitscore": ...}
        }
      }
    """
    cache_dir = Path(workdir) / "diamond_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    diamond_cache = {}

    for source_species in sorted(set(source_species_list)):
        if source_species == target_species:
            continue

        src_fa, tgt_fa = prepare_diamond_inputs(
            cache_dir=str(cache_dir),
            source_species=source_species,
            target_species=target_species,
            sequences_by_species=sequences_by_species,
        )

        out_tsv = cache_dir / f"{source_species}_vs_{target_species}.tsv"
        tmp_prefix = str(cache_dir / f"{source_species}_vs_{target_species}")

        if not out_tsv.exists():
            run_diamond(
                query_fa=src_fa,
                target_fa=tgt_fa,
                out_tsv=str(out_tsv),
                tmp_prefix=tmp_prefix,
            )

        diamond_cache[(source_species, target_species)] = load_diamond_results(str(out_tsv))

    return diamond_cache

def deduplicate_candidate_pairs(
    candidate_pairs: list[tuple[str, str, str, str, str, Interval | None]],
) -> list[tuple[str, str, str, str, str, Interval | None]]:
    seen = set()
    out = []

    for row in candidate_pairs:
        source_species, source_locus_id, target_species, target_locus_id, edge_origin, projected_interval = row

        proj_key = None
        if projected_interval is not None:
            proj_key = (
                projected_interval.seqid,
                projected_interval.start,
                projected_interval.end,
                projected_interval.strand,
            )

        key = (
            source_species,
            source_locus_id,
            target_species,
            target_locus_id,
            edge_origin,
            proj_key,
        )

        if key in seen:
            continue
        seen.add(key)
        out.append(row)

    return out

def deduplicate_candidate_pairs(
    candidate_pairs: list[tuple[str, str, str, str, str, Interval | None]],
) -> list[tuple[str, str, str, str, str, Interval | None]]:
    seen = set()
    out = []

    for row in candidate_pairs:
        source_species, source_locus_id, target_species, target_locus_id, edge_origin, projected_interval = row

        proj_key = None
        if projected_interval is not None:
            proj_key = (
                projected_interval.seqid,
                projected_interval.start,
                projected_interval.end,
                projected_interval.strand,
            )

        key = (
            source_species,
            source_locus_id,
            target_species,
            target_locus_id,
            edge_origin,
            proj_key,
        )

        if key in seen:
            continue
        seen.add(key)
        out.append(row)

    return out
    
def build_target_edge_evidence(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    hal_path: str,
    species_csv: str,
    merged_target_path: str,
) -> str:
    """
    Build orthology edge evidence for one target merged.json file.
    Returns the path to edge_evidence.json.
    """
    species_list = [x.strip() for x in species_csv.split(",") if x.strip()]
    transcripts_by_species = load_all_transcripts(annotation_dir, annotation_suffix, species_list)
    species_loci = build_all_species_loci(transcripts_by_species)

    merged_target = read_json(merged_target_path)
    target_species = merged_target["target_species"]

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

    candidate_pairs = deduplicate_candidate_pairs(candidate_pairs)
    
    seq_cache_dir = Path(workdir) / "sequence_cache"
    sequences_by_species = load_all_species_sequences(
        hal_path=hal_path,
        annotation_dir=annotation_dir,
        annotation_suffix=annotation_suffix,
        cache_dir=str(seq_cache_dir),
        species_list=species_list,
    )

    source_species_list = sorted({row[0] for row in candidate_pairs})
    diamond_cache = build_diamond_cache_for_target(
        workdir=workdir,
        target_species=target_species,
        source_species_list=source_species_list,
        sequences_by_species=sequences_by_species,
    )

    anchor_map = build_empty_anchor_map()

    edges, orthogroups = build_and_classify_edges(
        loci_by_species=species_loci,
        transcripts_by_species=build_transcript_lookup(transcripts_by_species),
        candidate_pairs=candidate_pairs,
        anchor_map=anchor_map,
        sequences_by_species=sequences_by_species,
        diamond_cache=diamond_cache,
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
