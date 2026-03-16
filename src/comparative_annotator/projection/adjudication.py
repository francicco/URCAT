from __future__ import annotations

from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.scoring import ProjectionLocusScore
from comparative_annotator.projection.scoring import score_projected_transcript_against_locus


def rank_locus_candidates(
    projected: ProjectedTranscript,
    loci: list[SpeciesLocus],
) -> list[ProjectionLocusScore]:

    scores = []

    for locus in loci:
        score = score_projected_transcript_against_locus(projected, locus)
        scores.append(score)

    scores.sort(key=lambda x: x.total_score, reverse=True)

    return scores


def choose_best_locus(
    projected: ProjectedTranscript,
    loci: list[SpeciesLocus],
):

    ranked = rank_locus_candidates(projected, loci)

    if not ranked:
        return None, []

    best = ranked[0]
    alternatives = ranked[1:]

    return best, alternatives
