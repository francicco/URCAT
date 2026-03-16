from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ReciprocalProjectionResult:
    source_species: str
    target_species: str

    seed_transcript: str
    forward_primary_locus: str | None
    reverse_primary_locus: str | None

    reciprocal_supported: bool
    classification: str
