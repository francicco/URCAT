from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ComparativeLocus:

    locus_id: str
    seed_species: str
    seed_transcript: str

    ortholog_loci: dict[str, list[str]] = field(default_factory=dict)
    paralog_loci: dict[str, list[str]] = field(default_factory=dict)

    missing_annotations: dict[str, list[str]] = field(default_factory=dict)

    def add_ortholog(self, species: str, locus_id: str):

        self.ortholog_loci.setdefault(species, []).append(locus_id)

    def add_paralog(self, species: str, locus_id: str):

        self.paralog_loci.setdefault(species, []).append(locus_id)

    def add_missing_projection(self, species: str, projection_id: str):

        self.missing_annotations.setdefault(species, []).append(projection_id)

    @property
    def species_count(self):

        return len(self.ortholog_loci)
