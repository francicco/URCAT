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

    # new fields
    source_exon_indices: list[int] = field(default_factory=list)
    chain_orientation: str | None = None   # "forward" or "reverse"
    fragmentation_count: int = 0
    chain_score: float | None = None

    def add_exon(self, start: int, end: int, source_exon_index: int | None = None) -> None:
        self.exons.append((start, end))
        if source_exon_index is not None:
            self.source_exon_indices.append(source_exon_index)

    @property
    def start(self) -> int:
        return min(e[0] for e in self.exons)

    @property
    def end(self) -> int:
        return max(e[1] for e in self.exons)

    @property
    def exon_count(self) -> int:
        return len(self.exons)
