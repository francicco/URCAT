from __future__ import annotations

from collections import defaultdict
from pathlib import Path


def _read_noncomment_lines(path: str | Path) -> list[str]:
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            out.append(line)
    return out


def _parse_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if not attr_str or attr_str == ".":
        return attrs

    for field in attr_str.split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value
        else:
            attrs[field] = ""
    return attrs


def _parse_gff_line(line: str, rec_index: int) -> dict:
    parts = line.split("\t")
    if len(parts) != 9:
        raise ValueError(f"Invalid GFF3 line with {len(parts)} columns: {line}")

    attrs = _parse_attributes(parts[8])

    return {
        "line": line,
        "index": rec_index,
        "seqid": parts[0],
        "source": parts[1],
        "type": parts[2],
        "start": int(parts[3]),
        "end": int(parts[4]),
        "score": parts[5],
        "strand": parts[6],
        "phase": parts[7],
        "attributes": attrs,
        "id": attrs.get("ID"),
        "parent": attrs.get("Parent"),
    }


def _feature_type_rank(feature_type: str) -> int:
    return {
        "gene": 0,
        "mRNA": 1,
        "transcript": 1,
        "exon": 2,
        "CDS": 3,
    }.get(feature_type, 9)


def _fallback_sort_key(rec: dict):
    return (
        rec["seqid"],
        rec["start"],
        rec["end"],
        _feature_type_rank(rec["type"]),
        rec["line"],
    )


def _transcript_child_sort_key(rec: dict):
    """
    Sort exon/CDS/features in transcript order.

    On + strand: ascending genomic coordinates.
    On - strand: descending genomic coordinates.
    For identical coordinates: exon before CDS.
    """
    strand = rec.get("strand", ".")
    start = rec["start"]
    end = rec["end"]
    type_rank = {"exon": 0, "CDS": 1}.get(rec["type"], 2)

    if strand == "-":
        return (-start, -end, type_rank, rec["index"])
    return (start, end, type_rank, rec["index"])


def _record_group_sort_key(rec: dict):
    return (
        rec["seqid"],
        rec["start"],
        rec["end"],
        _feature_type_rank(rec["type"]),
        rec["id"] or "",
        rec["index"],
    )


def _transcript_sort_key(tx_rec: dict, children: list[dict]):
    if children:
        child_starts = [c["start"] for c in children]
        child_ends = [c["end"] for c in children]
        tx_start = min([tx_rec["start"], *child_starts])
        tx_end = max([tx_rec["end"], *child_ends])
    else:
        tx_start = tx_rec["start"]
        tx_end = tx_rec["end"]

    strand = tx_rec.get("strand", ".")
    if strand == "-":
        return (tx_rec["seqid"], -tx_end, -tx_start, tx_rec["id"] or "", tx_rec["index"])
    return (tx_rec["seqid"], tx_start, tx_end, tx_rec["id"] or "", tx_rec["index"])


def _emit_transcript_block(tx_rec: dict, child_recs: list[dict]) -> list[str]:
    """
    Emit:
      transcript
      exon
      CDS
      exon
      CDS
      ...

    Any non-exon/CDS children are emitted after the transcript and before the
    exon/CDS stream.
    """
    out = [tx_rec["line"]]

    exon_recs = [r for r in child_recs if r["type"] == "exon"]
    cds_recs = [r for r in child_recs if r["type"] == "CDS"]
    other_recs = [r for r in child_recs if r["type"] not in {"exon", "CDS"}]

    other_recs = sorted(other_recs, key=_transcript_child_sort_key)
    out.extend(r["line"] for r in other_recs)

    exon_recs = sorted(exon_recs, key=_transcript_child_sort_key)
    cds_recs = sorted(cds_recs, key=_transcript_child_sort_key)

    cds_by_exact_coords: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for cds in cds_recs:
        cds_by_exact_coords[(cds["start"], cds["end"])].append(cds)

    used_cds_indices: set[int] = set()

    for exon in exon_recs:
        out.append(exon["line"])

        exact = cds_by_exact_coords.get((exon["start"], exon["end"]), [])
        exact = [c for c in exact if c["index"] not in used_cds_indices]
        if exact:
            for cds in exact:
                out.append(cds["line"])
                used_cds_indices.add(cds["index"])
            continue

        nested = [
            cds for cds in cds_recs
            if cds["index"] not in used_cds_indices
            and cds["start"] >= exon["start"]
            and cds["end"] <= exon["end"]
        ]
        nested = sorted(nested, key=_transcript_child_sort_key)
        for cds in nested:
            out.append(cds["line"])
            used_cds_indices.add(cds["index"])

    remaining_cds = [cds for cds in cds_recs if cds["index"] not in used_cds_indices]
    remaining_cds = sorted(remaining_cds, key=_transcript_child_sort_key)
    out.extend(r["line"] for r in remaining_cds)

    return out


def _order_gff_lines(lines: list[str]) -> list[str]:
    records = [_parse_gff_line(line, i) for i, line in enumerate(lines)]

    by_id: dict[str, dict] = {}
    for rec in records:
        if rec["id"]:
            by_id[rec["id"]] = rec

    gene_records: list[dict] = []
    transcript_records: list[dict] = []
    other_records: list[dict] = []

    for rec in records:
        if rec["type"] == "gene":
            gene_records.append(rec)
        elif rec["type"] in {"mRNA", "transcript"}:
            transcript_records.append(rec)
        else:
            other_records.append(rec)

    transcripts_by_gene: dict[str, list[dict]] = defaultdict(list)
    orphan_transcripts: list[dict] = []

    for tx in transcript_records:
        parent = tx.get("parent")
        if parent and parent in by_id and by_id[parent]["type"] == "gene":
            transcripts_by_gene[parent].append(tx)
        else:
            orphan_transcripts.append(tx)

    children_by_transcript: dict[str, list[dict]] = defaultdict(list)
    orphan_other_records: list[dict] = []

    for rec in other_records:
        parent = rec.get("parent")
        if parent and parent in by_id and by_id[parent]["type"] in {"mRNA", "transcript"}:
            children_by_transcript[parent].append(rec)
        else:
            orphan_other_records.append(rec)

    ordered_lines: list[str] = []

    gene_records = sorted(gene_records, key=_record_group_sort_key)
    for gene in gene_records:
        ordered_lines.append(gene["line"])

        txs = transcripts_by_gene.get(gene["id"], [])
        txs = sorted(
            txs,
            key=lambda tx: _transcript_sort_key(tx, children_by_transcript.get(tx["id"], [])),
        )

        for tx in txs:
            ordered_lines.extend(
                _emit_transcript_block(tx, children_by_transcript.get(tx["id"], []))
            )

    orphan_transcripts = sorted(
        orphan_transcripts,
        key=lambda tx: _transcript_sort_key(tx, children_by_transcript.get(tx["id"], [])),
    )
    for tx in orphan_transcripts:
        ordered_lines.extend(_emit_transcript_block(tx, children_by_transcript.get(tx["id"], [])))

    orphan_other_records = sorted(orphan_other_records, key=_fallback_sort_key)
    ordered_lines.extend(rec["line"] for rec in orphan_other_records)

    return ordered_lines


def write_final_species_gff3s(
    output_dir: str,
    annotation_paths: dict[str, str],
    species_list: list[str],
) -> None:
    """
    Merge original annotations with any URCAT new_loci GFF3s found under
    output_dir and write one final GFF3 per species.

    Output order:
      gene
      transcript/mRNA
      exon
      CDS
      exon
      CDS
      ...

    Exon/CDS order follows transcript order, including minus-strand transcripts.
    """
    out_dir = Path(output_dir) / "final_gff3"
    out_dir.mkdir(parents=True, exist_ok=True)

    for species in species_list:
        original_path_str = annotation_paths.get(species, "")
        original_path = Path(original_path_str) if original_path_str else None

        species_new_files = sorted(Path(output_dir).rglob(f"{species}.new_loci.gff3"))

        original_lines = (
            _read_noncomment_lines(original_path)
            if original_path is not None and original_path.exists()
            else []
        )

        new_lines: list[str] = []
        for p in species_new_files:
            new_lines.extend(_read_noncomment_lines(p))

        merged_lines = _order_gff_lines(original_lines + new_lines)
        ordered_new_lines = _order_gff_lines(new_lines)

        final_path = out_dir / f"{species}.final.gff3"
        new_only_path = out_dir / f"{species}.new_loci.gff3"

        with open(final_path, "w") as fh:
            fh.write("##gff-version 3\n")
            for line in merged_lines:
                fh.write(line + "\n")

        with open(new_only_path, "w") as fh:
            fh.write("##gff-version 3\n")
            for line in ordered_new_lines:
                fh.write(line + "\n")
