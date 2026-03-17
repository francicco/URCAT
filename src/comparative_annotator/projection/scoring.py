from dataclasses import dataclass
from typing import Optional


@dataclass
class ProjectionScore:
    species: str
    projected_id: str

    exon_recovery: float
    chain_completeness: float
    compactness: float
    orientation_consistency: float
    fragmentation_penalty: float

    projection_score: float


@dataclass
class ProjectionLocusScore(ProjectionScore):
    locus_id: Optional[str]
    overlap_fraction: float

    total_score: float


# -------------------------
# Feature computations
# -------------------------

def exon_recovery(projected, source_transcript) -> float:
    if not source_transcript.exons:
        return 0.0
    recovered = len(projected.source_exon_indices or [])
    return recovered / len(source_transcript.exons)


def chain_completeness(projected, source_transcript) -> float:
    idx = sorted(projected.source_exon_indices or [])
    if len(idx) <= 1:
        return 0.0

    adjacent = sum(1 for i in range(len(idx) - 1) if idx[i + 1] == idx[i] + 1)
    max_adjacent = len(source_transcript.exons) - 1
    if max_adjacent <= 0:
        return 0.0

    return adjacent / max_adjacent


def compactness(projected, source_transcript) -> float:
    if not projected.exons or not source_transcript.exons:
        return 0.0

    proj_start = min(s for s, e in projected.exons)
    proj_end = max(e for s, e in projected.exons)
    proj_span = proj_end - proj_start + 1

    src_start = min(s for s, e in source_transcript.exons)
    src_end = max(e for s, e in source_transcript.exons)
    src_span = src_end - src_start + 1

    if proj_span == 0 or src_span == 0:
        return 0.0

    return min(proj_span, src_span) / max(proj_span, src_span)


def orientation_consistency(projected) -> float:
    # do NOT penalize strongly
    return 1.0 if projected.chain_orientation == "forward" else 0.7


def fragmentation_penalty(projected, source_transcript) -> float:
    if not source_transcript.exons:
        return 0.0
    return min(1.0, projected.fragmentation_count / len(source_transcript.exons))


def locus_overlap_fraction(projected, locus) -> float:
    if not projected.exons:
        return 0.0

    proj_start = min(s for s, e in projected.exons)
    proj_end = max(e for s, e in projected.exons)

    overlap = max(0, min(proj_end, locus.end) - max(proj_start, locus.start))
    proj_len = proj_end - proj_start

    if proj_len <= 0:
        return 0.0

    return overlap / proj_len


# -------------------------
# Scoring functions
# -------------------------

def score_projected_transcript(projected, source_transcript) -> ProjectionScore:
    er = exon_recovery(projected, source_transcript)
    cc = chain_completeness(projected, source_transcript)
    comp = compactness(projected, source_transcript)
    ori = orientation_consistency(projected)
    frag = fragmentation_penalty(projected, source_transcript)

    projection_score = (
        0.4 * er +
        0.2 * cc +
        0.15 * comp +
        0.15 * ori -
        0.1 * frag
    )

    return ProjectionScore(
        species=projected.species,
        projected_id=projected.source_transcript,
        exon_recovery=er,
        chain_completeness=cc,
        compactness=comp,
        orientation_consistency=ori,
        fragmentation_penalty=frag,
        projection_score=projection_score,
    )


def score_projected_transcript_against_locus(projected, source_transcript, locus) -> ProjectionLocusScore:
    base = score_projected_transcript(projected, source_transcript)
    overlap = locus_overlap_fraction(projected, locus)

    total_score = (
        0.7 * base.projection_score +
        0.3 * overlap
    )

    return ProjectionLocusScore(
        **base.__dict__,
        locus_id=locus.locus_id if locus else None,
        overlap_fraction=overlap,
        total_score=total_score,
    )
