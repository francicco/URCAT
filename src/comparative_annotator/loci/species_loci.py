from __future__ import annotations

from collections import Counter
from typing import Iterable

from comparative_annotator.models.locus import SpeciesLocus


def infer_gene_id(transcript_ids):
    """
    Example:
        transcript:Eisa.Eisa2100G1.1 -> Eisa.Eisa2100G1
        transcript:Hmel202001oG3.1   -> Hmel202001oG3
    """
    cleaned = []
    for tx in transcript_ids:
        tx = tx.replace("transcript:", "")
        if "." in tx:
            tx = tx.rsplit(".", 1)[0]
        cleaned.append(tx)

    if not cleaned:
        return None

    counts = Counter(cleaned)
    return counts.most_common(1)[0][0]


def build_species_loci(transcripts: Iterable, species: str):
    """
    Group transcripts into species loci by genomic overlap on the same seqid/strand.
    Locus IDs are derived from the underlying transcript/gene IDs when possible,
    instead of synthetic labels like Eisa_locus_1.
    """
    txs = sorted(
        list(transcripts),
        key=lambda t: (t.seqid, t.strand, min(e[0] for e in t.exons), max(e[1] for e in t.exons)),
    )

    loci = []
    current = []

    def flush_current():
        nonlocal current, loci
        if not current:
            return

        seqid = current[0].seqid
        strand = current[0].strand
        start = min(min(e[0] for e in tx.exons) for tx in current)
        end = max(max(e[1] for e in tx.exons) for tx in current)
        transcript_ids = [tx.transcript_id for tx in current]

        gene_id = infer_gene_id(transcript_ids)
        if gene_id is None:
            gene_id = f"{species}_locus_{len(loci) + 1}"

        loci.append(
            SpeciesLocus(
                locus_id=gene_id,
                species=species,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                transcripts=transcript_ids,
            )
        )
        current = []

    for tx in txs:
        tx_start = min(e[0] for e in tx.exons)
        tx_end = max(e[1] for e in tx.exons)

        if not current:
            current = [tx]
            continue

        cur_seqid = current[0].seqid
        cur_strand = current[0].strand
        cur_end = max(max(e[1] for e in t.exons) for t in current)

        same_block = (tx.seqid == cur_seqid) and (tx.strand == cur_strand) and (tx_start <= cur_end)

        if same_block:
            current.append(tx)
        else:
            flush_current()
            current = [tx]

    flush_current()
    return loci
