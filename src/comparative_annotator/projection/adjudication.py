from __future__ import annotations

from comparative_annotator.projection.scoring import score_projected_transcript_against_locus


def rank_locus_candidates(
    projected,
    source_transcript,
    loci,
):
    scores = []

    for locus in loci:
        score = score_projected_transcript_against_locus(
            projected=projected,
            source_transcript=source_transcript,
            locus=locus,
        )
        scores.append(score)

    scores.sort(key=lambda x: x.total_score, reverse=True)
    return scores


def choose_best_locus(
    projected,
    source_transcript,
    loci,
):
    ranked = rank_locus_candidates(
        projected=projected,
        source_transcript=source_transcript,
        loci=loci,
    )

    if not ranked:
        return None, []

    return ranked[0], ranked[1:]
