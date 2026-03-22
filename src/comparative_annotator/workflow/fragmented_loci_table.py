from __future__ import annotations

import csv
from pathlib import Path

from comparative_annotator.workflow.fragmented_models import LogicalProjectedLocus


def write_fragmented_loci_table(
    out_path: str | Path,
    round_id: int,
    reference_species: str,
    species: str,
    loci: list[LogicalProjectedLocus],
) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cols = [
        "round_id",
        "reference_species",
        "species",
        "source_species",
        "locus_class",
        "locus_status",
        "n_target_seqids",
        "target_seqids",
        "dominant_seqid",
        "dominant_bp_fraction",
        "support_count",
        "total_chain_score",
        "mean_exon_recovery",
        "source_transcripts",
        "fragment_coords",
    ]

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()

        for locus in loci:
            if not locus.is_fragmented_across_seqids:
                continue

            fragment_coords = ",".join(
                f"{f.seqid}:{f.start}-{f.end}:{f.strand}" for f in locus.fragments
            )

            w.writerow(
                {
                    "round_id": round_id,
                    "reference_species": reference_species,
                    "species": species,
                    "source_species": locus.source_species,
                    "locus_class": locus.locus_class,
                    "locus_status": locus.locus_status,
                    "n_target_seqids": locus.n_target_seqids,
                    "target_seqids": ",".join(locus.seqids),
                    "dominant_seqid": locus.dominant_seqid,
                    "dominant_bp_fraction": f"{locus.dominant_bp_fraction:.3f}",
                    "support_count": locus.support_count,
                    "total_chain_score": f"{locus.total_chain_score:.3f}",
                    "mean_exon_recovery": f"{locus.mean_exon_recovery:.3f}",
                    "source_transcripts": ",".join(locus.source_transcripts),
                    "fragment_coords": fragment_coords,
                }
            )
