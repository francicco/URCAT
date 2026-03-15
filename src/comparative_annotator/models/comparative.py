from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ComparativeLocus:

    locus_id: str
    seed_species: str
    seed_transcript: str

    species_loci: dict[str, list[str]] = field(default_factory=dict)

    def add_species_locus(self, species: str, locus_id: str):

        if species not in self.species_loci:
            self.species_loci[species] = []

        if locus_id not in self.species_loci[species]:
            self.species_loci[species].append(locus_id)

    @property
    def species_count(self):

        return len(self.species_loci)
