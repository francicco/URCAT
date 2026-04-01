from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ComparativeLocus:
    """
    Canonical comparative locus container for URCAT.

    This object stores, for one seed transcript, the inferred locus-level outcome
    in each target species:
      - primary locus match
      - alternative locus matches
      - missing annotation / no convincing annotated locus
      - strand conflicts

    Notes
    -----
    - `primary` and `alternatives` store locus IDs.
    - `primary_transcripts` and `alternative_transcripts` store transcript IDs.
    - `missing_annotations` stores projection IDs or source transcript IDs that
      were not explained by an annotated target locus.
    - `strand_conflicts` stores locus IDs that overlapped but on the wrong strand.
    """

    locus_id: str
    seed_species: str
    seed_transcript: str

    primary: dict[str, str] = field(default_factory=dict)
    primary_transcripts: dict[str, str] = field(default_factory=dict)

    alternatives: dict[str, list[str]] = field(default_factory=dict)
    alternative_transcripts: dict[str, list[str]] = field(default_factory=dict)

    missing_annotations: dict[str, list[str]] = field(default_factory=dict)
    strand_conflicts: dict[str, list[str]] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Primary calls
    # ------------------------------------------------------------------

    def set_primary(self, species: str, locus_id: str) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        if not locus_id:
            raise ValueError("locus_id must be a non-empty string")
        self.primary[species] = locus_id

    def set_primary_transcript(self, species: str, transcript_id: str) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        if not transcript_id:
            raise ValueError("transcript_id must be a non-empty string")
        self.primary_transcripts[species] = transcript_id

    # ------------------------------------------------------------------
    # Alternative calls
    # ------------------------------------------------------------------

    def set_alternatives(self, species: str, locus_ids: list[str]) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        self.alternatives[species] = self._dedup_preserve_order(locus_ids)

    def add_alternative(self, species: str, locus_id: str) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        if not locus_id:
            raise ValueError("locus_id must be a non-empty string")
        self.alternatives.setdefault(species, [])
        if locus_id not in self.alternatives[species]:
            self.alternatives[species].append(locus_id)

    def set_alternative_transcripts(self, species: str, transcript_ids: list[str]) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        self.alternative_transcripts[species] = self._dedup_preserve_order(transcript_ids)

    def add_alternative_transcript(self, species: str, transcript_id: str) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        if not transcript_id:
            raise ValueError("transcript_id must be a non-empty string")
        self.alternative_transcripts.setdefault(species, [])
        if transcript_id not in self.alternative_transcripts[species]:
            self.alternative_transcripts[species].append(transcript_id)

    # ------------------------------------------------------------------
    # Missing / conflict annotations
    # ------------------------------------------------------------------

    def add_missing_projection(self, species: str, projection_id: str) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        if not projection_id:
            raise ValueError("projection_id must be a non-empty string")
        self.missing_annotations.setdefault(species, [])
        if projection_id not in self.missing_annotations[species]:
            self.missing_annotations[species].append(projection_id)

    def add_strand_conflict(self, species: str, locus_id: str) -> None:
        if not species:
            raise ValueError("species must be a non-empty string")
        if not locus_id:
            raise ValueError("locus_id must be a non-empty string")
        self.strand_conflicts.setdefault(species, [])
        if locus_id not in self.strand_conflicts[species]:
            self.strand_conflicts[species].append(locus_id)

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------

    def has_primary(self, species: str) -> bool:
        return species in self.primary and bool(self.primary[species])

    def has_alternatives(self, species: str) -> bool:
        return species in self.alternatives and len(self.alternatives[species]) > 0

    def has_missing(self, species: str) -> bool:
        return species in self.missing_annotations and len(self.missing_annotations[species]) > 0

    def has_strand_conflict(self, species: str) -> bool:
        return species in self.strand_conflicts and len(self.strand_conflicts[species]) > 0

    @property
    def species_count(self) -> int:
        """
        Number of species with a primary locus assignment.
        """
        return len(self.primary)

    def to_dict(self) -> dict:
        return {
            "locus_id": self.locus_id,
            "seed_species": self.seed_species,
            "seed_transcript": self.seed_transcript,
            "primary": self.primary,
            "primary_transcripts": self.primary_transcripts,
            "alternatives": self.alternatives,
            "alternative_transcripts": self.alternative_transcripts,
            "missing_annotations": self.missing_annotations,
            "strand_conflicts": self.strand_conflicts,
        }

    @staticmethod
    def _dedup_preserve_order(values: list[str] | None) -> list[str]:
        if not values:
            return []
        seen = set()
        out = []
        for value in values:
            if not value:
                continue
            if value not in seen:
                seen.add(value)
                out.append(value)
        return out
