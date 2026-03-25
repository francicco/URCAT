from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class FragmentedBlock:
    source_exon_number: int
    target_seqid: str
    target_start: int
    target_end: int
    target_strand: str
    chain_score: float | None = None


@dataclass
class FragmentedComparativeLocus:
    source_species: str
    source_transcript: str
    target_species: str

    source_seqid: str
    source_strand: str
    n_source_exons: int

    target_seqids: list[str]
    blocks: list[FragmentedBlock] = field(default_factory=list)

    recovered_exon_numbers: list[int] = field(default_factory=list)
    exon_recovery_fraction: float = 0.0
    strand_consistent: bool = True
    exon_order_consistent: bool = True
    total_chain_score: float = 0.0

    join_score: float = 0.0
    status: str = "fragmented"

    aa_identity_mean: float | None = None
    aa_coverage_fraction: float | None = None
    aa_bitscore: float | None = None
    protein_support_class: str | None = None

    fragment_class: str | None = None
    classification_reason: str | None = None

    max_consecutive_recovered_exons: int | None = None
    consecutive_exon_fraction: float | None = None

    def to_row(self) -> dict:
        return {
            "source_species": self.source_species,
            "source_transcript": self.source_transcript,
            "target_species": self.target_species,
            "source_seqid": self.source_seqid,
            "source_strand": self.source_strand,
            "n_source_exons": self.n_source_exons,
            "target_seqids": ",".join(self.target_seqids),
            "n_target_seqids": len(self.target_seqids),
            "recovered_exon_numbers": ",".join(map(str, self.recovered_exon_numbers)),
            "n_recovered_exons": len(self.recovered_exon_numbers),
            "exon_recovery_fraction": round(self.exon_recovery_fraction, 4),
            "strand_consistent": self.strand_consistent,
            "exon_order_consistent": self.exon_order_consistent,
            "total_chain_score": round(self.total_chain_score, 4),
            "join_score": round(self.join_score, 4),
            "aa_identity_mean": None if self.aa_identity_mean is None else round(self.aa_identity_mean, 4),
            "aa_coverage_fraction": None if self.aa_coverage_fraction is None else round(self.aa_coverage_fraction, 4),
            "aa_bitscore": None if self.aa_bitscore is None else round(self.aa_bitscore, 4),
            "protein_support_class": self.protein_support_class,
            "fragment_class": self.fragment_class,
            "classification_reason": self.classification_reason,
            "status": self.status,
        }
