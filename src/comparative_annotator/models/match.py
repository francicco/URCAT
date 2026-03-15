from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ProjectionLocusMatch:
    projected_id: str
    species: str
    locus_id: str
    overlap_fraction: float
    classification: str
