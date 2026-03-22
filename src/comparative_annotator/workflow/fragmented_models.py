from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ProjectedBlock:
    source_transcript_id: str
    source_exon_index: int
    target_species: str
    target_seqid: str
    target_start: int
    target_end: int
    target_strand: str
    chain_score: float | None = None
    exon_recovery: float | None = None


@dataclass
class FragmentModel:
    seqid: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(x[0] for x in self.exons)

    @property
    def end(self) -> int:
        return max(x[1] for x in self.exons)

    @property
    def span(self) -> int:
        return self.end - self.start + 1

    @property
    def exon_bp(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)


@dataclass
class LogicalProjectedLocus:
    species: str
    source_species: str | None
    source_transcripts: list[str]
    support_count: int
    total_chain_score: float
    mean_exon_recovery: float
    fragments: list[FragmentModel]
    dominant_seqid: str
    dominant_bp_fraction: float
    n_target_seqids: int
    is_fragmented_across_seqids: bool
    locus_class: str
    locus_status: str
    notes: list[str] = field(default_factory=list)

    @property
    def seqids(self) -> list[str]:
        return [f.seqid for f in self.fragments]

    @property
    def total_exonic_bp(self) -> int:
        return sum(f.exon_bp for f in self.fragments)

    def to_dict(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "source_species": self.source_species,
            "source_transcripts": self.source_transcripts,
            "support_count": self.support_count,
            "total_chain_score": self.total_chain_score,
            "mean_exon_recovery": self.mean_exon_recovery,
            "dominant_seqid": self.dominant_seqid,
            "dominant_bp_fraction": self.dominant_bp_fraction,
            "n_target_seqids": self.n_target_seqids,
            "is_fragmented_across_seqids": self.is_fragmented_across_seqids,
            "locus_class": self.locus_class,
            "locus_status": self.locus_status,
            "notes": self.notes,
            "fragments": [
                {
                    "seqid": f.seqid,
                    "strand": f.strand,
                    "start": f.start,
                    "end": f.end,
                    "exons": f.exons,
                    "exon_bp": f.exon_bp,
                }
                for f in self.fragments
            ],
        }
