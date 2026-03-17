from __future__ import annotations


def write_merged_gff3(
    output_path,
    original_transcripts_by_species,
    species,
    urcat_consensus_loci,
):
    """
    Write a merged GFF3:
    - original annotation
    - URCAT-inferred missing loci appended at the end

    Parameters
    ----------
    output_path : str
    original_transcripts_by_species : dict[str, dict[str, Transcript]]
    species : str
    urcat_consensus_loci : list[MissingLocusConsensus]
    """

    with open(output_path, "w") as out:
        # Header
        out.write("##gff-version 3\n")

        # -------------------------
        # 1. ORIGINAL ANNOTATION
        # -------------------------
        transcripts = original_transcripts_by_species.get(species, {})

        for tx in transcripts.values():
            gene_id = tx.transcript_id.split(".")[0]  # crude but OK for now

            # gene
            out.write(
                "\t".join([
                    tx.seqid,
                    "annotation",
                    "gene",
                    str(tx.start),
                    str(tx.end),
                    ".",
                    tx.strand,
                    ".",
                    f"ID={gene_id}",
                ]) + "\n"
            )

            # mRNA
            out.write(
                "\t".join([
                    tx.seqid,
                    "annotation",
                    "mRNA",
                    str(tx.start),
                    str(tx.end),
                    ".",
                    tx.strand,
                    ".",
                    f"ID={tx.transcript_id};Parent={gene_id}",
                ]) + "\n"
            )

            # exons
            for i, (start, end) in enumerate(tx.exons, start=1):
                exon_id = f"{tx.transcript_id}.exon{i}"
                out.write(
                    "\t".join([
                        tx.seqid,
                        "annotation",
                        "exon",
                        str(start),
                        str(end),
                        ".",
                        tx.strand,
                        ".",
                        f"ID={exon_id};Parent={tx.transcript_id}",
                    ]) + "\n"
                )

        # -------------------------
        # 2. URCAT NEW LOCI
        # -------------------------
        counter = 1

        for consensus in urcat_consensus_loci:
            gene_id = f"URCAT_{species}_new_{counter:06d}"
            tx_id = f"{gene_id}.t1"

            start = min(s for s, e in consensus.exons)
            end = max(e for s, e in consensus.exons)

            attributes_common = (
                f"urcat_status=new_locus;"
                f"urcat_support_count={consensus.support_count};"
                f"urcat_total_chain_score={consensus.total_chain_score:.3f};"
                f"urcat_source_transcripts={','.join(consensus.source_transcripts)};"
                f"Note=Projected missing locus inferred by URCAT"
            )

            # gene
            out.write(
                "\t".join([
                    consensus.seqid,
                    "URCAT",
                    "gene",
                    str(start),
                    str(end),
                    ".",
                    consensus.strand,
                    ".",
                    f"ID={gene_id};{attributes_common}",
                ]) + "\n"
            )

            # mRNA
            out.write(
                "\t".join([
                    consensus.seqid,
                    "URCAT",
                    "mRNA",
                    str(start),
                    str(end),
                    ".",
                    consensus.strand,
                    ".",
                    f"ID={tx_id};Parent={gene_id};{attributes_common}",
                ]) + "\n"
            )

            # exons
            for i, (start, end) in enumerate(consensus.exons, start=1):
                exon_id = f"{tx_id}.exon{i}"
                out.write(
                    "\t".join([
                        consensus.seqid,
                        "URCAT",
                        "exon",
                        str(start),
                        str(end),
                        ".",
                        consensus.strand,
                        ".",
                        f"ID={exon_id};Parent={tx_id}",
                    ]) + "\n"
                )

            counter += 1
