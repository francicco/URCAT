from __future__ import annotations

from comparative_annotator.projection.scoring import score_projected_transcript_against_locus
from comparative_annotator.projection.transcript_ranking import rank_transcripts_within_locus
from comparative_annotator.projection.matching import find_candidate_species_loci

from comparative_annotator.projection.matching import (
    find_candidate_species_loci,
    locus_relation_to_projection,
    locus_distance_to_projection,
)

def rank_candidate_loci_with_transcripts(
    projected,
    source_transcript,
    species_loci,
    target_transcripts_by_id=None,
):
    if target_transcripts_by_id is None:
        target_transcripts_by_id = {}

    candidate_loci = find_candidate_species_loci(
        projected,
        species_loci,
        n_flank=2,
    )

    rows = []

    for locus in candidate_loci:
        locus_score = score_projected_transcript_against_locus(
            projected=projected,
            source_transcript=source_transcript,
            locus=locus,
        )

        transcript_scores = rank_transcripts_within_locus(
            projected=projected,
            source_transcript=source_transcript,
            locus=locus,
            target_transcripts_by_id=target_transcripts_by_id,
        )

        best_tx = transcript_scores[0] if transcript_scores else None

        rows.append(
            {
                "locus_id": locus.locus_id,
                "locus_strand": locus.strand,
                "projected_strand": projected.strand,
                "same_strand": locus.strand == projected.strand,
                "locus_score": locus_score.total_score,
                "overlap_fraction": locus_score.overlap_fraction,
                "projection_score": locus_score.projection_score,
                "best_transcript_id": None if best_tx is None else best_tx.transcript_id,
                "best_transcript_score": None if best_tx is None else best_tx.total_score,
                "transcript_scores": transcript_scores,
            }
        )

    rows.sort(
        key=lambda x: (
            x["locus_score"],
            -int(not x["same_strand"]),
            x["best_transcript_score"] if x["best_transcript_score"] is not None else -1.0,
        ),
        reverse=True,
    )

    return rows
