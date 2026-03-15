from comparative_annotator.models.locus import SpeciesLocus


def build_species_loci(transcripts):

    txs = sorted(
        transcripts,
        key=lambda x: (x.seqid, x.start)
    )

    loci = []

    current = None
    counter = 0

    for tx in txs:

        if current is None:

            counter += 1

            current = SpeciesLocus(
                locus_id=f"locus_{counter}",
                species=tx.species,
                seqid=tx.seqid,
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                transcripts=[tx.transcript_id],
            )

            continue

        if tx.start <= current.end:

            current.transcripts.append(tx.transcript_id)
            current.end = max(current.end, tx.end)

        else:

            loci.append(current)

            counter += 1

            current = SpeciesLocus(
                locus_id=f"locus_{counter}",
                species=tx.species,
                seqid=tx.seqid,
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                transcripts=[tx.transcript_id],
            )

    if current:
        loci.append(current)

    return loci
