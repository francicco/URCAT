from __future__ import annotations


def _attrs_to_str(attrs: dict[str, str]) -> str:
    parts = []
    for key, value in attrs.items():
        parts.append(f"{key}={value}")
    return ";".join(parts)


def _write_feature(
    out,
    seqid,
    source,
    feature_type,
    start,
    end,
    strand,
    attributes,
    score=".",
    phase=".",
):
    out.write(
        "\t".join(
            [
                str(seqid),
                str(source),
                str(feature_type),
                str(start),
                str(end),
                str(score),
                str(strand),
                str(phase),
                _attrs_to_str(attributes),
            ]
        )
        + "\n"
    )


def _guess_gene_id_from_transcript_id(transcript_id: str) -> str:
    """
    Very simple fallback.
    You may later want to replace this with a parser using real Parent/gene IDs.
    """
    if ".t" in transcript_id:
        return transcript_id.split(".t")[0]
    if transcript_id.startswith("transcript:"):
        return transcript_id.replace("transcript:", "gene:", 1)
    return f"{transcript_id}.gene"


def write_full_merged_gff3(
    output_path,
    original_transcripts_by_species,
    urcat_consensus_by_species=None,
):
    """
    Write one merged GFF3 containing:
    - original annotation for all species
    - appended URCAT-inferred loci for any species with new loci

    Parameters
    ----------
    output_path : str
    original_transcripts_by_species : dict[str, dict[str, CandidateTranscript]]
        Example:
            {
                "Hmel": {...},
                "Eisa": {...},
                "Diul": {...},
            }

    urcat_consensus_by_species : dict[str, list[MissingLocusConsensus]] | None
        Example:
            {
                "Eisa": [consensus1, consensus2],
                "Diul": [],
            }
    """
    if urcat_consensus_by_species is None:
        urcat_consensus_by_species = {}

    with open(output_path, "w") as out:
        out.write("##gff-version 3\n")

        # -------------------------
        # 1. ORIGINAL ANNOTATION
        # -------------------------
        for species in sorted(original_transcripts_by_species):
            transcripts = original_transcripts_by_species[species]

            out.write(f"### ORIGINAL ANNOTATION: {species}\n")

            # Stable order
            for tx_id in sorted(transcripts):
                tx = transcripts[tx_id]

                gene_id = _guess_gene_id_from_transcript_id(tx.transcript_id)

                # gene
                _write_feature(
                    out=out,
                    seqid=tx.seqid,
                    source="annotation",
                    feature_type="gene",
                    start=tx.start,
                    end=tx.end,
                    strand=tx.strand,
                    attributes={
                        "ID": gene_id,
                        "species": species,
                    },
                )

                # mRNA
                _write_feature(
                    out=out,
                    seqid=tx.seqid,
                    source="annotation",
                    feature_type="mRNA",
                    start=tx.start,
                    end=tx.end,
                    strand=tx.strand,
                    attributes={
                        "ID": tx.transcript_id,
                        "Parent": gene_id,
                        "species": species,
                    },
                )

                # exons
                for i, (start, end) in enumerate(tx.exons, start=1):
                    exon_id = f"{tx.transcript_id}.exon{i}"
                    _write_feature(
                        out=out,
                        seqid=tx.seqid,
                        source="annotation",
                        feature_type="exon",
                        start=start,
                        end=end,
                        strand=tx.strand,
                        attributes={
                            "ID": exon_id,
                            "Parent": tx.transcript_id,
                            "species": species,
                        },
                    )

            out.write("###\n")

        # -------------------------
        # 2. URCAT NEW LOCI
        # -------------------------
        for species in sorted(urcat_consensus_by_species):
            consensus_loci = urcat_consensus_by_species.get(species, [])
            if not consensus_loci:
                continue

            out.write(f"### URCAT NEW LOCI: {species}\n")

            for counter, consensus in enumerate(consensus_loci, start=1):
                gene_id = f"URCAT_{species}_new_{counter:06d}"
                tx_id = f"{gene_id}.t1"

                start = min(s for s, e in consensus.exons)
                end = max(e for s, e in consensus.exons)

                common_attrs = {
                    "species": species,
                    "urcat_status": "new_locus",
                    "urcat_support_count": str(consensus.support_count),
                    "urcat_total_chain_score": f"{consensus.total_chain_score:.3f}",
                    "urcat_mean_exon_recovery": f"{consensus.mean_exon_recovery:.3f}",
                    "urcat_source_transcripts": ",".join(consensus.source_transcripts),
                    "Note": "Projected missing locus inferred by URCAT",
                }

                # gene
                gene_attrs = {"ID": gene_id, **common_attrs}
                _write_feature(
                    out=out,
                    seqid=consensus.seqid,
                    source="URCAT",
                    feature_type="gene",
                    start=start,
                    end=end,
                    strand=consensus.strand,
                    attributes=gene_attrs,
                )

                # mRNA
                mrna_attrs = {"ID": tx_id, "Parent": gene_id, **common_attrs}
                _write_feature(
                    out=out,
                    seqid=consensus.seqid,
                    source="URCAT",
                    feature_type="mRNA",
                    start=start,
                    end=end,
                    strand=consensus.strand,
                    attributes=mrna_attrs,
                )

                # exons
                for i, (exon_start, exon_end) in enumerate(consensus.exons, start=1):
                    exon_id = f"{tx_id}.exon{i}"
                    exon_attrs = {
                        "ID": exon_id,
                        "Parent": tx_id,
                        "species": species,
                    }
                    _write_feature(
                        out=out,
                        seqid=consensus.seqid,
                        source="URCAT",
                        feature_type="exon",
                        start=exon_start,
                        end=exon_end,
                        strand=consensus.strand,
                        attributes=exon_attrs,
                    )

            out.write("###\n")
