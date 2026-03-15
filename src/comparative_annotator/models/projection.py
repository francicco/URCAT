from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ProjectionInterval:
    species: str
    seqid: str
    start: int
    end: int
    strand: str

    source_species: str
    source_transcript: str

    identity: float | None = None
    coverage: float | None = None
