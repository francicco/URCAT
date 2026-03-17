from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ProjectionTranscriptScore:
    species: str
    projected_id: str
    transcript_id: str
    locus_id: str | None

    exon_recovery: float
    exon_count_similarity: float
    span_overlap: float
    orientation_consistency: float
    intron_support: float

    total_score: float
