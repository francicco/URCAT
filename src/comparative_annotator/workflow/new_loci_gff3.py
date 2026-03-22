from __future__ import annotations

from pathlib import Path

from comparative_annotator.workflow.fragmented_models import LogicalProjectedLocus


def _attrs_to_str(attrs: dict[str, object]) -> str:
    out = []
    for k, v in attrs.items():
        if v is None:
            continue
        out.append(f"{k}={v}")
    return ";".join(out)


def _feature_line(
    seqid: str,
    source: str,
    feature_type: str,
    start: int,
    end: int,
    strand: str,
    attrs: dict[str, object],
    score: str = ".",
    phase: str = ".",
) -> str:
    return "\t".join(
        [
            seqid,
            source,
            feature_type,
            str(start),
            str(end),
            score,
            strand,
            phase,
            _attrs_to_str(attrs),
        ]
    )


def _logical_gene_base_id(species: str, idx: int) -> str:
    return f"{species}.URCAT.newG{idx:06d}"


def _fragment_group_id(species: str, idx: int) -> str:
    return f"{species}.URCAT.fragG{idx:06d}"


def write_new_loci_gff3(
    out_path: str | Path,
    species: str,
    loci: list[LogicalProjectedLocus],
) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")

        gene_index = 1
        for locus in sorted(
            loci,
            key=lambda x: (
                x.fragments[0].seqid if x.fragments else "",
                min(f.start for f in x.fragments) if x.fragments else 0,
                max(f.end for f in x.fragments) if x.fragments else 0,
            ),
        ):
            if locus.is_fragmented_across_seqids:
                frag_group = _fragment_group_id(species, gene_index)
                for frag_i, frag in enumerate(
                    sorted(locus.fragments, key=lambda f: (f.seqid, f.start, f.end)),
                    start=1,
                ):
                    gene_id = f"{frag_group}.f{frag_i}"
                    tx_id = f"{gene_id}.1"

                    common_attrs = {
                        "source": "URCAT",
                        "urcat_status": locus.locus_status,
                        "urcat_class": locus.locus_class,
                        "urcat_fragment_group": frag_group,
                        "urcat_fragment_index": f"{frag_i}/{len(locus.fragments)}",
                        "urcat_n_seqids": locus.n_target_seqids,
                        "urcat_target_seqids": ",".join(locus.seqids),
                        "urcat_dominant_seqid": locus.dominant_seqid,
                        "urcat_dominant_bp_fraction": f"{locus.dominant_bp_fraction:.3f}",
                        "urcat_support_count": locus.support_count,
                        "urcat_total_chain_score": f"{locus.total_chain_score:.3f}",
                        "urcat_mean_exon_recovery": f"{locus.mean_exon_recovery:.3f}",
                        "urcat_source_transcripts": ",".join(locus.source_transcripts),
                    }

                    fh.write(
                        _feature_line(
                            frag.seqid,
                            "URCAT",
                            "gene",
                            frag.start,
                            frag.end,
                            frag.strand,
                            {"ID": gene_id, "Name": gene_id, **common_attrs},
                        )
                        + "\n"
                    )
                    fh.write(
                        _feature_line(
                            frag.seqid,
                            "URCAT",
                            "mRNA",
                            frag.start,
                            frag.end,
                            frag.strand,
                            {"ID": tx_id, "Parent": gene_id, "Name": tx_id, **common_attrs},
                        )
                        + "\n"
                    )

                    for exon_i, (s, e) in enumerate(frag.exons, start=1):
                        fh.write(
                            _feature_line(
                                frag.seqid,
                                "URCAT",
                                "exon",
                                s,
                                e,
                                frag.strand,
                                {"ID": f"{tx_id}.exon{exon_i}", "Parent": tx_id},
                            )
                            + "\n"
                        )
                        fh.write(
                            _feature_line(
                                frag.seqid,
                                "URCAT",
                                "CDS",
                                s,
                                e,
                                frag.strand,
                                {"ID": f"{tx_id}.cds{exon_i}", "Parent": tx_id},
                                phase=".",
                            )
                            + "\n"
                        )
                gene_index += 1
                continue

            frag = locus.fragments[0]
            gene_id = _logical_gene_base_id(species, gene_index)
            tx_id = f"{gene_id}.1"

            common_attrs = {
                "source": "URCAT",
                "urcat_status": locus.locus_status,
                "urcat_class": locus.locus_class,
                "urcat_support_count": locus.support_count,
                "urcat_total_chain_score": f"{locus.total_chain_score:.3f}",
                "urcat_mean_exon_recovery": f"{locus.mean_exon_recovery:.3f}",
                "urcat_source_transcripts": ",".join(locus.source_transcripts),
            }

            fh.write(
                _feature_line(
                    frag.seqid,
                    "URCAT",
                    "gene",
                    frag.start,
                    frag.end,
                    frag.strand,
                    {"ID": gene_id, "Name": gene_id, **common_attrs},
                )
                + "\n"
            )
            fh.write(
                _feature_line(
                    frag.seqid,
                    "URCAT",
                    "mRNA",
                    frag.start,
                    frag.end,
                    frag.strand,
                    {"ID": tx_id, "Parent": gene_id, "Name": tx_id, **common_attrs},
                )
                + "\n"
            )

            for exon_i, (s, e) in enumerate(frag.exons, start=1):
                fh.write(
                    _feature_line(
                        frag.seqid,
                        "URCAT",
                        "exon",
                        s,
                        e,
                        frag.strand,
                        {"ID": f"{tx_id}.exon{exon_i}", "Parent": tx_id},
                    )
                    + "\n"
                )
                fh.write(
                    _feature_line(
                        frag.seqid,
                        "URCAT",
                        "CDS",
                        s,
                        e,
                        frag.strand,
                        {"ID": f"{tx_id}.cds{exon_i}", "Parent": tx_id},
                        phase=".",
                    )
                    + "\n"
                )

            gene_index += 1
