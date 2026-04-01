from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ComparativeLocus:
    seed_species: str
    seed_transcript: str

    primary: dict[str, str] = field(default_factory=dict)
    alternatives: dict[str, list[str]] = field(default_factory=dict)
    missing_annotations: dict[str, list[str]] = field(default_factory=dict)
    strand_conflicts: dict[str, list[str]] = field(default_factory=dict)

    primary_transcripts: dict[str, str] = field(default_factory=dict)
    alternative_transcripts: dict[str, list[str]] = field(default_factory=dict)

    def set_primary(self, species: str, locus_id: str, transcript_id: str | None = None) -> None:
        self.primary[species] = locus_id
        if transcript_id:
            self.primary_transcripts[species] = transcript_id

    def add_alternative(self, species: str, locus_id: str, transcript_id: str | None = None) -> None:
        self.alternatives.setdefault(species, [])
        if locus_id not in self.alternatives[species]:
            self.alternatives[species].append(locus_id)

        if transcript_id:
            self.alternative_transcripts.setdefault(species, [])
            if transcript_id not in self.alternative_transcripts[species]:
                self.alternative_transcripts[species].append(transcript_id)

    def add_missing_projection(self, species: str, projection_id: str) -> None:
        self.missing_annotations.setdefault(species, [])
        if projection_id not in self.missing_annotations[species]:
            self.missing_annotations[species].append(projection_id)

    def add_strand_conflict(self, species: str, projection_id: str) -> None:
        self.strand_conflicts.setdefault(species, [])
        if projection_id not in self.strand_conflicts[species]:
            self.strand_conflicts[species].append(projection_id)
