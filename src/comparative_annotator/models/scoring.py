from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ProjectionLocusScore:
    projected_id: str
    locus_id: str
    species: str

    overlap_fraction: float
    exon_count_similarity: float
    exon_recovery_fraction: float
    projection_coverage: float

    total_score: float
