from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ProjectedTranscript:

    species: str
    seqid: str
    strand: str

    source_species: str
    source_transcript: str

    exons: list[tuple[int, int]] = field(default_factory=list)

    identity: float | None = None
    coverage: float | None = None

    def add_exon(self, start: int, end: int):

        self.exons.append((start, end))

    @property
    def start(self):

        return min(e[0] for e in self.exons)

    @property
    def end(self):

        return max(e[1] for e in self.exons)

    @property
    def exon_count(self):

        return len(self.exons)
