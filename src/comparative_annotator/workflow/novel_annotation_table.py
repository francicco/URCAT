from __future__ import annotations

from pathlib import Path

from comparative_annotator.workflow.reporting import read_json, write_tsv


NOVEL_COLUMNS = [
    "round_id",
    "reference_species",
    "species",
    "event_type",
    "selected_for_next_round",
    "seqid",
    "start",
    "end",
    "strand",
    "support_count",
    "total_chain_score",
    "mean_exon_recovery",
    "source_transcripts",
    "source_loci",
    "locus_id",
    "transcripts",
]


def _interval_dict_to_row(interval: dict) -> dict:
    return {
        "seqid": interval.get("seqid"),
        "start": interval.get("start"),
        "end": interval.get("end"),
        "strand": interval.get("strand"),
    }


def build_novel_annotation_rows(round_dir: str | Path) -> list[dict]:
    round_dir = Path(round_dir)
    decision_path = round_dir / "post_round_decision.json"
    if not decision_path.exists():
        return []

    x = read_json(decision_path)
    round_id = x.get("round_id")
    reference_species = x.get("reference_species")
    next_reference_species = x.get("next_reference_species")

    rows = []

    for species, items in (x.get("new_consensus_by_species") or {}).items():
        for item in items:
            row = {
                "round_id": round_id,
                "reference_species": reference_species,
                "species": species,
                "event_type": "new_consensus",
                "selected_for_next_round": species == next_reference_species,
                "support_count": item.get("support_count"),
                "total_chain_score": item.get("total_chain_score"),
                "mean_exon_recovery": item.get("mean_exon_recovery"),
                "source_transcripts": item.get("source_transcripts"),
                "source_loci": item.get("source_loci"),
                "locus_id": item.get("locus_id"),
                "transcripts": item.get("transcripts"),
            }
            row.update(_interval_dict_to_row(item))
            rows.append(row)

    for species, items in (x.get("orphan_loci_by_species") or {}).items():
        for item in items:
            row = {
                "round_id": round_id,
                "reference_species": reference_species,
                "species": species,
                "event_type": "orphan_native_locus",
                "selected_for_next_round": species == next_reference_species,
                "support_count": None,
                "total_chain_score": None,
                "mean_exon_recovery": None,
                "source_transcripts": None,
                "source_loci": None,
                "locus_id": item.get("locus_id"),
                "transcripts": item.get("transcripts"),
            }
            row.update(_interval_dict_to_row(item))
            rows.append(row)

    for species, items in (x.get("pending_seeds_by_species") or {}).items():
        for item in items:
            row = {
                "round_id": round_id,
                "reference_species": reference_species,
                "species": species,
                "event_type": "pending_seed",
                "selected_for_next_round": species == next_reference_species,
                "support_count": item.get("support_count"),
                "total_chain_score": item.get("total_chain_score"),
                "mean_exon_recovery": item.get("mean_exon_recovery"),
                "source_transcripts": item.get("source_transcripts"),
                "source_loci": item.get("source_loci"),
                "locus_id": item.get("locus_id"),
                "transcripts": item.get("transcripts"),
            }
            row.update(_interval_dict_to_row(item))
            rows.append(row)

    return rows


def write_novel_annotations_table(round_dir: str | Path) -> str:
    round_dir = Path(round_dir)
    rows = build_novel_annotation_rows(round_dir)
    out_path = round_dir / "novel_annotations.tsv"
    write_tsv(out_path, rows, columns=NOVEL_COLUMNS)
    return str(out_path)
