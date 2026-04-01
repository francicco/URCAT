from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ProjectedTranscript:
    """
    Canonical projected transcript model for URCAT.

    Semantics
    ---------
    species
        Target species where this transcript was projected.
    source_species
        Species of the original seed transcript.
    source_transcript
        Identifier of the original source transcript.
    seqid
        Target seqid/scaffold/chromosome.
    strand
        Target strand.
    exons
        Target exon coordinates as 1-based inclusive intervals.

    Notes
    -----
    This must be the only ProjectedTranscript definition used across the codebase.
    Older code may expect fields such as exon_count, coverage, start, and end, so
    they are exposed as properties here.
    """

    species: str
    source_species: str
    source_transcript: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    chain_orientation: str | None = None
    chain_score: float | None = None
    mean_exon_recovery: float | None = None
    support_score: float | None = None

    def __post_init__(self) -> None:
        self.exons = [self._normalize_interval(start, end) for start, end in self.exons]
        self.exons.sort(key=lambda x: (x[0], x[1]))

    @staticmethod
    def _normalize_interval(start: int, end: int) -> tuple[int, int]:
        start_i = int(start)
        end_i = int(end)
        return (start_i, end_i) if start_i <= end_i else (end_i, start_i)

    @property
    def start(self) -> int:
        if not self.exons:
            raise ValueError("ProjectedTranscript has no exons, so start is undefined")
        return min(start for start, _ in self.exons)

    @property
    def end(self) -> int:
        if not self.exons:
            raise ValueError("ProjectedTranscript has no exons, so end is undefined")
        return max(end for _, end in self.exons)

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def coverage(self) -> int:
        """
        Spliced exon span in nucleotides, using 1-based inclusive coordinates.
        """
        return sum((end - start + 1) for start, end in self.exons)

    @property
    def intron_count(self) -> int:
        return max(0, len(self.exons) - 1)

    def span_tuple(self) -> tuple[str, int, int, str]:
        return (self.seqid, self.start, self.end, self.strand)

    def same_intron_chain(self, other: "ProjectedTranscript") -> bool:
        return (
            self.seqid == other.seqid
            and self.strand == other.strand
            and self.exons == other.exons
        )

    def overlaps(self, other: "ProjectedTranscript") -> bool:
        if self.seqid != other.seqid:
            return False
        return not (self.end < other.start or other.end < self.start)

    def to_dict(self) -> dict:
        return {
            "species": self.species,
            "source_species": self.source_species,
            "source_transcript": self.source_transcript,
            "seqid": self.seqid,
            "strand": self.strand,
            "start": self.start if self.exons else None,
            "end": self.end if self.exons else None,
            "exons": self.exons,
            "exon_count": self.exon_count,
            "coverage": self.coverage,
            "chain_orientation": self.chain_orientation,
            "chain_score": self.chain_score,
            "mean_exon_recovery": self.mean_exon_recovery,
            "support_score": self.support_score,
        }
