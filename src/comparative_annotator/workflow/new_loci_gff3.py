from __future__ import annotations

from pathlib import Path


def _attr(obj, name, default=None):
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def _consensus_sort_key(x):
    seqid = _attr(x, "seqid", "")
    exons = _attr(x, "exons", []) or []
    start = min((e[0] for e in exons), default=10**18)
    end = max((e[1] for e in exons), default=-1)
    strand = _attr(x, "strand", ".")
    return (seqid, start, end, strand)


def _format_attrs(attrs):
    return ";".join(f"{k}={v}" for k, v in attrs.items() if v is not None)


def write_new_loci_gff3(out_path: str, target_species: str, loci) -> str:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    gene_counter = 0

    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")

        for locus in sorted(loci, key=_consensus_sort_key):
            exons = list(_attr(locus, "exons", []) or [])
            if not exons:
                continue

            gene_counter += 1
            gene_id = f"{target_species}.URCAT.newG{gene_counter:06d}"
            tx_id = f"{gene_id}.1"

            seqid = _attr(locus, "seqid")
            strand = _attr(locus, "strand", ".")
            start = min(e[0] for e in exons)
            end = max(e[1] for e in exons)

            support_count = _attr(locus, "support_count")
            total_chain_score = _attr(locus, "total_chain_score")
            mean_exon_recovery = _attr(locus, "mean_exon_recovery")
            source_transcripts = _attr(locus, "source_transcripts", [])
            if isinstance(source_transcripts, (list, tuple)):
                source_transcripts = ",".join(source_transcripts)

            gene_attrs = {
                "ID": gene_id,
                "Name": gene_id,
                "source": "URCAT",
                "urcat_status": "new_locus",
                "urcat_support_count": (
                    f"{support_count}" if support_count is not None else None
                ),
                "urcat_total_chain_score": (
                    f"{total_chain_score:.3f}"
                    if isinstance(total_chain_score, (int, float))
                    else total_chain_score
                ),
                "urcat_mean_exon_recovery": (
                    f"{mean_exon_recovery:.3f}"
                    if isinstance(mean_exon_recovery, (int, float))
                    else mean_exon_recovery
                ),
                "urcat_source_transcripts": source_transcripts if source_transcripts else None,
            }

            tx_attrs = {
                "ID": tx_id,
                "Parent": gene_id,
                "Name": tx_id,
                "source": "URCAT",
                "urcat_status": "new_locus",
                "urcat_support_count": (
                    f"{support_count}" if support_count is not None else None
                ),
                "urcat_total_chain_score": (
                    f"{total_chain_score:.3f}"
                    if isinstance(total_chain_score, (int, float))
                    else total_chain_score
                ),
                "urcat_mean_exon_recovery": (
                    f"{mean_exon_recovery:.3f}"
                    if isinstance(mean_exon_recovery, (int, float))
                    else mean_exon_recovery
                ),
                "urcat_source_transcripts": source_transcripts if source_transcripts else None,
            }

            fh.write(
                "\t".join(
                    [
                        str(seqid),
                        "URCAT",
                        "gene",
                        str(start),
                        str(end),
                        ".",
                        str(strand),
                        ".",
                        _format_attrs(gene_attrs),
                    ]
                )
                + "\n"
            )

            fh.write(
                "\t".join(
                    [
                        str(seqid),
                        "URCAT",
                        "mRNA",
                        str(start),
                        str(end),
                        ".",
                        str(strand),
                        ".",
                        _format_attrs(tx_attrs),
                    ]
                )
                + "\n"
            )

            exon_order = sorted(exons, key=lambda e: (e[0], e[1]))
            if strand == "-":
                exon_order = list(reversed(exon_order))

            for i, (exon_start, exon_end) in enumerate(exon_order, start=1):
                exon_attrs = {
                    "ID": f"{tx_id}.exon{i}",
                    "Parent": tx_id,
                }
                fh.write(
                    "\t".join(
                        [
                            str(seqid),
                            "URCAT",
                            "exon",
                            str(exon_start),
                            str(exon_end),
                            ".",
                            str(strand),
                            ".",
                            _format_attrs(exon_attrs),
                        ]
                    )
                    + "\n"
                )

    return str(out_path)
