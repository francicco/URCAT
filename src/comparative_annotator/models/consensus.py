from dataclasses import dataclass


@dataclass
class ComparativeConsensus:

    expected_cds_length: float | None = None
    expected_exon_count: float | None = None
