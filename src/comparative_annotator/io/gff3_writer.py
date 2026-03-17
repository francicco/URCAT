from __future__ import annotations

from pathlib import Path


def _attrs_to_str(attrs: dict[str, str]) -> str:
    return ";".join(f"{k}={v}" for k, v in attrs.items())


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


def _get_tx_attr(tx, name, default=None):
    return getattr(tx, name, default)


def _original_gene_id(tx) -> str:
    gene_id = _get_tx_attr(tx, "gene_id", None)
    if gene_id:
        return gene_id

    attrs = _get_tx_attr(tx, "attributes", None)
    if isinstance(attrs, dict):
        if "gene_id" in attrs:
            return attrs["gene_id"]
        if "gene" in attrs:
            return attrs["gene"]
        if "Parent" in attrs and attrs["Parent"]:
            return attrs["Parent"]

    return f"{tx.transcript_id}.gene"


def _original_transcript_id(tx) -> str:
    tx_id = _get_tx_attr(tx, "transcript_id", None)
    if tx_id:
        return tx_id

    attrs = _get_tx_attr(tx, "attributes", None)
    if isinstance(attrs, dict) and "ID" in attrs:
        return attrs["ID"]

    raise ValueError("Transcript object has no transcript_id")


def write_species_merged_gff3(
    output_path,
    transcripts_by_species,
    species,
    urcat_consensus_loci=None,
):
    """
    Write one GFF3 for a single species:
    - original annotation with preserved IDs
    - optional appended URCAT-inferred loci
    """
    if urcat_consensus_loci is None:
        urcat_consensus_loci = []

    transcripts = transcripts_by_species.get(species, {})

    with open(output_path, "w") as out:
        out.write("##gff-version 3\n")
        out.write(f"### ORIGINAL ANNOTATION: {species}\n")

        for tx_id in sorted(transcripts):
            tx = transcripts[tx_id]

            gene_id = _original_gene_id(tx)
            transcript_id = _original_transcript_id(tx)
            source = _get_tx_attr(tx, "source", "annotation")

            # gene
            _write_feature(
                out=out,
                seqid=tx.seqid,
                source=source,
                feature_type="gene",
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                attributes={
                    "ID": gene_id,
                },
            )

            # mRNA
            _write_feature(
                out=out,
                seqid=tx.seqid,
                source=source,
                feature_type="mRNA",
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                attributes={
                    "ID": transcript_id,
                    "Parent": gene_id,
                },
            )

            # exons
            for i, (start, end) in enumerate(tx.exons, start=1):
                exon_id = f"{transcript_id}.exon{i}"
                _write_feature(
                    out=out,
                    seqid=tx.seqid,
                    source=source,
                    feature_type="exon",
                    start=start,
                    end=end,
                    strand=tx.strand,
                    attributes={
                        "ID": exon_id,
                        "Parent": transcript_id,
                    },
                )

        out.write("###\n")

        if urcat_consensus_loci:
            out.write(f"### URCAT NEW LOCI: {species}\n")

        for counter, consensus in enumerate(urcat_consensus_loci, start=1):
            gene_id = f"URCAT_{species}_new_{counter:06d}"
            transcript_id = f"{gene_id}.t1"

            start = min(s for s, e in consensus.exons)
            end = max(e for s, e in consensus.exons)

            common_attrs = {
                "urcat_status": "new_locus",
                "urcat_support_count": str(consensus.support_count),
                "urcat_total_chain_score": f"{consensus.total_chain_score:.3f}",
                "urcat_mean_exon_recovery": f"{consensus.mean_exon_recovery:.3f}",
                "urcat_source_transcripts": ",".join(consensus.source_transcripts),
                "Note": "Projected missing locus inferred by URCAT",
            }

            _write_feature(
                out=out,
                seqid=consensus.seqid,
                source="URCAT",
                feature_type="gene",
                start=start,
                end=end,
                strand=consensus.strand,
                attributes={
                    "ID": gene_id,
                    **common_attrs,
                },
            )

            _write_feature(
                out=out,
                seqid=consensus.seqid,
                source="URCAT",
                feature_type="mRNA",
                start=start,
                end=end,
                strand=consensus.strand,
                attributes={
                    "ID": transcript_id,
                    "Parent": gene_id,
                    **common_attrs,
                },
            )

            for i, (exon_start, exon_end) in enumerate(consensus.exons, start=1):
                exon_id = f"{transcript_id}.exon{i}"
                _write_feature(
                    out=out,
                    seqid=consensus.seqid,
                    source="URCAT",
                    feature_type="exon",
                    start=exon_start,
                    end=exon_end,
                    strand=consensus.strand,
                    attributes={
                        "ID": exon_id,
                        "Parent": transcript_id,
                    },
                )

        if urcat_consensus_loci:
            out.write("###\n")


def write_species_gff3_outputs(
    output_dir,
    transcripts_by_species,
    urcat_consensus_by_species=None,
):
    """
    Write one merged GFF3 per species.
    Returns dict: species -> output path
    """
    if urcat_consensus_by_species is None:
        urcat_consensus_by_species = {}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}

    for species in sorted(transcripts_by_species):
        outpath = output_dir / f"{species}.merged.urcat.gff3"
        write_species_merged_gff3(
            output_path=str(outpath),
            transcripts_by_species=transcripts_by_species,
            species=species,
            urcat_consensus_loci=urcat_consensus_by_species.get(species, []),
        )
        outputs[species] = str(outpath)

    return outputs
