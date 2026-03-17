from __future__ import annotations

from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.transcript import CandidateTranscript


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def build_species_loci(
    transcripts: list[CandidateTranscript],
    species: str,
) -> list[SpeciesLocus]:
    """
    Build species-local loci by grouping overlapping transcripts
    on the same seqid and strand.
    """

    def infer_gene_id(transcript_ids):
        # example transcript: "transcript:Eisa.Eisa2100G1.1"
        # -> gene: "Eisa.Eisa2100G1"
        tx = transcript_ids[0]
        tx = tx.replace("transcript:", "")
        return tx.rsplit(".", 1)[0]

    txs = sorted(
        [tx for tx in transcripts if tx.species == species],
        key=lambda x: (x.seqid, x.strand, x.start, x.end),
    )

    loci: list[SpeciesLocus] = []
    current: SpeciesLocus | None = None
    counter = 0

    for tx in txs:
        if current is None:
            counter += 1
            current = SpeciesLocus(
                locus_id = infer_gene_id(locus.transcripts)
                species=species,
                seqid=tx.seqid,
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                transcripts=[tx.transcript_id],
            )
            tx.attributes["locus_id"] = current.locus_id
            continue

        same_context = (
            tx.seqid == current.seqid
            and tx.strand == current.strand
        )

        if same_context and overlaps(tx.start, tx.end, current.start, current.end):
            current.transcripts.append(tx.transcript_id)
            current.start = min(current.start, tx.start)
            current.end = max(current.end, tx.end)
            tx.attributes["locus_id"] = current.locus_id
        else:
            loci.append(current)
            counter += 1
            current = SpeciesLocus(
                locus_id=f"{species}_locus_{counter}",
                species=species,
                seqid=tx.seqid,
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                transcripts=[tx.transcript_id],
            )
            tx.attributes["locus_id"] = current.locus_id

    if current is not None:
        loci.append(current)

    return loci
