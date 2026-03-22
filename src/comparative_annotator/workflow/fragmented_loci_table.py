from __future__ import annotations

import csv
from pathlib import Path


def _get(row, key, default=""):
    value = row.get(key, default)
    if value is None:
        return default
    if isinstance(value, list):
        return ",".join(str(x) for x in value)
    return value


def write_fragmented_loci_table(
    out_path: str | Path,
    round_id: int,
    reference_species: str,
    species: str,
    loci: list[dict],
) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "round_id",
        "reference_species",
        "species",
        "source_species",
        "source_transcripts",
        "target_species",
        "is_fragmented",
        "is_fragmented_across_seqids",
        "n_target_seqids",
        "target_seqids",
        "n_blocks",
        "block_sizes",
        "target_span_start",
        "target_span_end",
        "target_span_bp",
    ]

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for locus in loci:
            is_fragmented = bool(
                locus.get("is_fragmented")
                or locus.get("is_fragmented_across_seqids")
                or int(locus.get("n_target_seqids", 1)) > 1
            )

            if not is_fragmented:
                continue

            writer.writerow(
                {
                    "round_id": round_id,
                    "reference_species": reference_species,
                    "species": species,
                    "source_species": _get(locus, "source_species"),
                    "source_transcripts": _get(locus, "source_transcripts"),
                    "target_species": _get(locus, "target_species", species),
                    "is_fragmented": _get(locus, "is_fragmented", True),
                    "is_fragmented_across_seqids": _get(
                        locus, "is_fragmented_across_seqids", False
                    ),
                    "n_target_seqids": _get(locus, "n_target_seqids"),
                    "target_seqids": _get(locus, "target_seqids"),
                    "n_blocks": _get(locus, "n_blocks"),
                    "block_sizes": _get(locus, "block_sizes"),
                    "target_span_start": _get(locus, "target_span_start"),
                    "target_span_end": _get(locus, "target_span_end"),
                    "target_span_bp": _get(locus, "target_span_bp"),
                }
            )
