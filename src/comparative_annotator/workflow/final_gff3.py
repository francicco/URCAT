from __future__ import annotations

from pathlib import Path


def _parse_attrs(text: str) -> dict[str, str]:
    out = {}
    for field in text.strip().split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        k, v = field.split("=", 1)
        out[k] = v
    return out


def _format_attrs(attrs: dict[str, str]) -> str:
    priority = ["ID", "Parent", "Name"]
    items = []

    for k in priority:
        if k in attrs:
            items.append((k, attrs[k]))
    for k in sorted(attrs):
        if k not in {"ID", "Parent", "Name"}:
            items.append((k, attrs[k]))

    return ";".join(f"{k}={v}" for k, v in items)


def _read_gff3(path: str | Path) -> list[dict]:
    path = Path(path)
    rows = []

    if not path.exists():
        return rows

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            rows.append(
                {
                    "seqid": seqid,
                    "source": source,
                    "type": ftype,
                    "start": int(start),
                    "end": int(end),
                    "score": score,
                    "strand": strand,
                    "phase": phase,
                    "attrs": _parse_attrs(attrs),
                }
            )

    return rows


def _feature_rank(ftype: str) -> int:
    return {
        "gene": 0,
        "mRNA": 1,
        "transcript": 1,
        "exon": 2,
        "CDS": 3,
        "five_prime_UTR": 4,
        "three_prime_UTR": 5,
    }.get(ftype, 99)


def _deduplicate_rows(rows: list[dict]) -> list[dict]:
    seen = set()
    out = []

    for r in rows:
        key = (
            r["seqid"],
            r["source"],
            r["type"],
            r["start"],
            r["end"],
            r["score"],
            r["strand"],
            r["phase"],
            tuple(sorted(r["attrs"].items())),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(r)

    return out


def _sort_rows(rows: list[dict]) -> list[dict]:
    return sorted(
        rows,
        key=lambda r: (
            r["seqid"],
            r["start"],
            r["end"],
            _feature_rank(r["type"]),
            r["attrs"].get("ID", ""),
            r["attrs"].get("Parent", ""),
        ),
    )


def _write_gff3(path: str | Path, rows: list[dict]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    rows = _sort_rows(_deduplicate_rows(rows))

    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for r in rows:
            out.write(
                "\t".join(
                    [
                        r["seqid"],
                        r["source"],
                        r["type"],
                        str(r["start"]),
                        str(r["end"]),
                        r["score"],
                        r["strand"],
                        r["phase"],
                        _format_attrs(r["attrs"]),
                    ]
                )
                + "\n"
            )
    return str(path)


def collect_round_new_loci_gff3(output_dir: str | Path, species: str) -> list[str]:
    output_dir = Path(output_dir)
    hits = []

    for p in sorted(output_dir.rglob(f"{species}.new_loci.gff3")):
        if p.exists() and p.stat().st_size > 0:
            hits.append(str(p))

    return hits


def write_species_new_loci_gff3(output_dir: str | Path, species: str) -> str:
    output_dir = Path(output_dir)
    out_dir = output_dir / "final_gff3"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in collect_round_new_loci_gff3(output_dir, species):
        rows.extend(_read_gff3(p))

    out_path = out_dir / f"{species}.new_loci.gff3"
    return _write_gff3(out_path, rows)


def write_species_final_gff3(
    output_dir: str | Path,
    annotation_dir: str | Path,
    annotation_suffix: str,
    species: str,
) -> str:
    output_dir = Path(output_dir)
    annotation_dir = Path(annotation_dir)
    out_dir = output_dir / "final_gff3"
    out_dir.mkdir(parents=True, exist_ok=True)

    native_path = annotation_dir / f"{species}{annotation_suffix}"
    new_loci_path = out_dir / f"{species}.new_loci.gff3"

    rows = []
    rows.extend(_read_gff3(native_path))
    rows.extend(_read_gff3(new_loci_path))

    out_path = out_dir / f"{species}.final.gff3"
    return _write_gff3(out_path, rows)


def write_all_final_species_gff3(
    output_dir: str | Path,
    annotation_dir: str | Path,
    annotation_suffix: str,
    species_list: list[str],
) -> list[str]:
    written = []
    for species in species_list:
        write_species_new_loci_gff3(output_dir, species)
        written.append(
            write_species_final_gff3(
                output_dir=output_dir,
                annotation_dir=annotation_dir,
                annotation_suffix=annotation_suffix,
                species=species,
            )
        )
    return written
