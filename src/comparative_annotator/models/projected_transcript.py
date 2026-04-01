from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ProjectedTranscript:
    """
    Canonical projected-transcript model used across URCAT.

    Semantics
    ---------
    species:
        Target species where the projection lands.
    source_species:
        Source/reference species from which the transcript was projected.
    source_transcript:
        Source transcript identifier.
    seqid / strand / exons:
        Target-side reconstructed transcript structure.
    """

    species: str
    source_species: str
    source_transcript: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    target_seqids: list[str] = field(default_factory=list)
    chain_orientation: str | None = None
    total_chain_score: float = 0.0
    mean_chain_score: float = 0.0
    exon_recovery_fraction: float = 0.0
    n_source_exons: int = 0
    recovered_exon_numbers: list[int] = field(default_factory=list)

    source_seqid: str | None = None
    source_strand: str | None = None

    @property
    def start(self) -> int:
        return min(s for s, _ in self.exons) if self.exons else 0

    @property
    def end(self) -> int:
        return max(e for _, e in self.exons) if self.exons else 0

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def coverage(self) -> float:
        return self.exon_recovery_fraction

    def to_dict(self) -> dict:
        return {
            "species": self.species,
            "source_species": self.source_species,
            "source_transcript": self.source_transcript,
            "seqid": self.seqid,
            "strand": self.strand,
            "exons": [list(x) for x in self.exons],
            "target_seqids": list(self.target_seqids),
            "chain_orientation": self.chain_orientation,
            "total_chain_score": self.total_chain_score,
            "mean_chain_score": self.mean_chain_score,
            "exon_recovery_fraction": self.exon_recovery_fraction,
            "n_source_exons": self.n_source_exons,
            "recovered_exon_numbers": list(self.recovered_exon_numbers),
            "source_seqid": self.source_seqid,
            "source_strand": self.source_strand,
        }
