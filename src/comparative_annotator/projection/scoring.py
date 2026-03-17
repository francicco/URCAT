from __future__ import annotations

from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.scoring import ProjectionLocusScore
from comparative_annotator.projection.matching import (
    locus_overlap_fraction,
    projected_transcript_id,
)


def exon_count_similarity(projected: ProjectedTranscript, locus: SpeciesLocus) -> float:
    if projected.exon_count == 0:
        return 0.0
    return 1.0


def exon_recovery_fraction(projected: ProjectedTranscript) -> float:
    if projected.coverage is None:
        return 0.0
    # normalize roughly by projected exon count
    return float(projected.coverage) / max(1, projected.exon_count)


def projection_coverage_score(projected: ProjectedTranscript) -> float:
    if projected.chain_score is not None:
        # compress chain score into something modest
        return max(0.0, min(1.0, projected.chain_score / max(1, projected.exon_count)))
    if projected.coverage is None:
        return 0.0
    return 1.0


def score_projected_transcript_against_locus(
    projected: ProjectedTranscript,
    locus: SpeciesLocus,
    w_overlap: float = 0.5,
    w_exon_similarity: float = 0.15,
    w_exon_recovery: float = 0.2,
    w_projection_coverage: float = 0.15,
) -> ProjectionLocusScore:
    overlap = locus_overlap_fraction(projected, locus)
    exon_sim = exon_count_similarity(projected, locus)
    recovery = exon_recovery_fraction(projected)
    proj_cov = projection_coverage_score(projected)

    total = (
        w_overlap * overlap
        + w_exon_similarity * exon_sim
        + w_exon_recovery * recovery
        + w_projection_coverage * proj_cov
    )

    return ProjectionLocusScore(
        projected_id=projected_transcript_id(projected),
        locus_id=locus.locus_id,
        species=projected.species,
        overlap_fraction=overlap,
        exon_count_similarity=exon_sim,
        exon_recovery_fraction=recovery,
        projection_coverage=proj_cov,
        total_score=total,
    )
