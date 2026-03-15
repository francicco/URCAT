from dataclasses import dataclass, field


@dataclass
class SpeciesLocus:

    locus_id: str
    species: str
    seqid: str
    start: int
    end: int
    strand: str

    transcripts: list[str] = field(default_factory=list)
