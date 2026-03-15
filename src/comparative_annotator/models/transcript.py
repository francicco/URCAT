from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class CandidateTranscript:
    transcript_id: str
    species: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str

    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)

    cds_length: int = 0
    exon_count: int = 0
    intron_count: int = 0
    cdna_length: int = 0

    local_score: float | None = None
    projection_score: float | None = None
    comparative_score: float | None = None
    total_score: float | None = None

    attributes: dict[str, Any] = field(default_factory=dict)

    def finalize(self) -> None:
        self.exons.sort()
        self.cds.sort()

        self.exon_count = len(self.exons)
        self.intron_count = max(0, self.exon_count - 1)
        self.cdna_length = sum(e2 - e1 + 1 for e1, e2 in self.exons)
        self.cds_length = sum(c2 - c1 + 1 for c1, c2 in self.cds)

        if self.exons:
            self.start = min(e1 for e1, _ in self.exons)
            self.end = max(e2 for _, e2 in self.exons)
        elif self.cds:
            self.start = min(c1 for c1, _ in self.cds)
            self.end = max(c2 for _, c2 in self.cds)

    @property
    def intron_chain(self) -> list[tuple[int, int]]:
        if len(self.exons) < 2:
            return []

        introns = []
        for i in range(len(self.exons) - 1):
            intron_start = self.exons[i][1] + 1
            intron_end = self.exons[i + 1][0] - 1
            introns.append((intron_start, intron_end))

        return introns
