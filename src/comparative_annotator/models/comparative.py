from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ComparativeLocus:
    locus_id: str
    seed_species: str
    seed_transcript: str

    primary: dict[str, str] = field(default_factory=dict)
    alternatives: dict[str, list[str]] = field(default_factory=dict)
    missing_annotations: dict[str, list[str]] = field(default_factory=dict)
    strand_conflicts: dict[str, list[str]] = field(default_factory=dict)

    def set_primary(self, species: str, locus_id: str) -> None:
        self.primary[species] = locus_id

    def set_alternatives(self, species: str, locus_ids: list[str]) -> None:
        self.alternatives[species] = locus_ids

    def add_missing_projection(self, species: str, projection_id: str) -> None:
        self.missing_annotations.setdefault(species, []).append(projection_id)

    def add_strand_conflict(self, species: str, locus_id: str) -> None:
        self.strand_conflicts.setdefault(species, []).append(locus_id)

    @property
    def species_count(self) -> int:
        return len(self.primary)
