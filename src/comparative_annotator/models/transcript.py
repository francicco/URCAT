from dataclasses import dataclass, field


@dataclass
class CandidateTranscript:

    transcript_id: str
    species: str
    seqid: str
    start: int
    end: int
    strand: str

    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)

    cds_length: int = 0
    exon_count: int = 0

    local_score: float | None = None
    projection_score: float | None = None
    comparative_score: float | None = None
    total_score: float | None = None

    def finalize(self):

        self.exons.sort()
        self.cds.sort()

        self.exon_count = len(self.exons)

        self.cds_length = sum(
            e2 - e1 + 1 for e1, e2 in self.cds
        )
