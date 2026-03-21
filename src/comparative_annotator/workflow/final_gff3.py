from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class NewLocusModel:
    species: str
    seqid: str
    start: int
    end: int
    strand: str
    exons: list[tuple[int, int]]
    support_count: int | None = None
    total_chain_score: float | None = None
    mean_exon_recovery: float | None = None
    source_transcripts: list[str] | None = None


def _format_attrs(attrs: dict[str, Any]) -> str:
    parts = []
    for k, v in attrs.items():
        if v is None or v == "":
            continue
        parts.append(f"{k}={v}")
    return ";".join(parts)


def _locus_sort_key(x: NewLocusModel):
    return (
        x.seqid,
        x.start,
        x.end,
        x.strand,
        ",".join(x.source_transcripts or []),
    )


def _feature_rank(feature_type: str) -> int:
    return {
        "gene": 0,
        "mRNA": 1,
        "transcript": 1,
        "exon": 2,
        "CDS": 3,
    }.get(feature_type, 9)


def _safe_int(x: Any, default: int = 0) -> int:
    try:
        return int(x)
    except Exception:
        return default


def _safe_float(x: Any, default: float | None = None) -> float | None:
    try:
        return float(x)
    except Exception:
        return default


def _normalize_exons(exons: list[Any]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for e in exons or []:
        if not isinstance(e, (list, tuple)) or len(e) != 2:
            continue
        s = _safe_int(e[0], -1)
        t = _safe_int(e[1], -1)
        if s <= 0 or t <= 0:
            continue
        if t < s:
            s, t = t, s
        out.append((s, t))
    out.sort()
    return out


def _deduplicate_new_loci(loci: list[NewLocusModel]) -> list[NewLocusModel]:
    seen = set()
    out = []

    for loc in sorted(loci, key=_locus_sort_key):
        key = (
            loc.species,
            loc.seqid,
            loc.start,
            loc.end,
            loc.strand,
            tuple(loc.exons),
            tuple(loc.source_transcripts or []),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(loc)

    return out


def renumber_new_loci(species: str, loci: list[NewLocusModel]) -> list[dict[str, Any]]:
    """
    Assign globally unique URCAT IDs for one species across all rounds.
    """
    loci = _deduplicate_new_loci(loci)
    loci_sorted = sorted(loci, key=_locus_sort_key)

    out = []
    for i, loc in enumerate(loci_sorted, start=1):
        gene_id = f"{species}.URCAT.newG{i:06d}"
        tx_id = f"{gene_id}.1"
        exon_ids = [f"{tx_id}.exon{j}" for j in range(1, len(loc.exons) + 1)]

        out.append(
            {
                "model": loc,
                "gene_id": gene_id,
                "tx_id": tx_id,
                "exon_ids": exon_ids,
            }
        )

    return out


def render_new_loci_gff3_lines(species: str, loci: list[NewLocusModel]) -> list[str]:
    """
    Render renumbered URCAT loci as GFF3 lines, keeping each locus block together.
    """
    lines: list[str] = []
    renumbered = renumber_new_loci(species, loci)

    for item in renumbered:
        loc = item["model"]
        gene_id = item["gene_id"]
        tx_id = item["tx_id"]
        exon_ids = item["exon_ids"]

        common_attrs = {
            "source": "URCAT",
            "urcat_status": "new_locus",
            "urcat_support_count": loc.support_count,
            "urcat_total_chain_score": (
                f"{loc.total_chain_score:.3f}" if loc.total_chain_score is not None else None
            ),
            "urcat_mean_exon_recovery": (
                f"{loc.mean_exon_recovery:.3f}" if loc.mean_exon_recovery is not None else None
            ),
            "urcat_source_transcripts": ",".join(loc.source_transcripts or []),
        }

        gene_attrs = {
            "ID": gene_id,
            "Name": gene_id,
            **common_attrs,
        }
        lines.append(
            "\t".join(
                [
                    loc.seqid,
                    "URCAT",
                    "gene",
                    str(loc.start),
                    str(loc.end),
                    ".",
                    loc.strand,
                    ".",
                    _format_attrs(gene_attrs),
                ]
            )
        )

        tx_attrs = {
            "ID": tx_id,
            "Parent": gene_id,
            "Name": tx_id,
            **common_attrs,
        }
        lines.append(
            "\t".join(
                [
                    loc.seqid,
                    "URCAT",
                    "mRNA",
                    str(loc.start),
                    str(loc.end),
                    ".",
                    loc.strand,
                    ".",
                    _format_attrs(tx_attrs),
                ]
            )
        )

        for j, (ex_start, ex_end) in enumerate(loc.exons):
            exon_attrs = {
                "ID": exon_ids[j],
                "Parent": tx_id,
            }
            lines.append(
                "\t".join(
                    [
                        loc.seqid,
                        "URCAT",
                        "exon",
                        str(ex_start),
                        str(ex_end),
                        ".",
                        loc.strand,
                        ".",
                        _format_attrs(exon_attrs),
                    ]
                )
            )

    return lines


def write_new_loci_gff3(path: str | Path, species: str, loci: list[NewLocusModel]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for line in render_new_loci_gff3_lines(species, loci):
            out.write(line + "\n")


def parse_gff3_records(gff_path: str | Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    with open(gff_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = cols
            start_i = _safe_int(start, 0)
            end_i = _safe_int(end, 0)

            records.append(
                {
                    "seqid": seqid,
                    "start": start_i,
                    "end": end_i,
                    "feature": feature,
                    "strand": strand,
                    "line": line,
                }
            )

    return records


def merge_original_and_new_gff3(
    original_gff3: str | Path,
    final_gff3: str | Path,
    species: str,
    new_loci: list[NewLocusModel],
) -> None:
    final_gff3 = Path(final_gff3)
    final_gff3.parent.mkdir(parents=True, exist_ok=True)

    original_records = parse_gff3_records(original_gff3)
    new_lines = render_new_loci_gff3_lines(species, new_loci)

    new_records = []
    for line in new_lines:
        cols = line.split("\t")
        seqid, source, feature, start, end, score, strand, phase, attrs = cols
        new_records.append(
            {
                "seqid": seqid,
                "start": _safe_int(start, 0),
                "end": _safe_int(end, 0),
                "feature": feature,
                "strand": strand,
                "line": line,
            }
        )

    all_records = original_records + new_records
    all_records.sort(
        key=lambda r: (
            r["seqid"],
            r["start"],
            r["end"],
            _feature_rank(r["feature"]),
            r["strand"],
            r["line"],
        )
    )

    with open(final_gff3, "w") as out:
        out.write("##gff-version 3\n")
        for rec in all_records:
            out.write(rec["line"] + "\n")


def load_all_new_loci_for_species(rounds_dir: str | Path, species: str) -> list[NewLocusModel]:
    """
    Pool all new consensuses for one species across all rounds.
    """
    rounds_dir = Path(rounds_dir)
    loci: list[NewLocusModel] = []

    for decision_path in sorted(rounds_dir.glob("round_*/post_round_decision.json")):
        with open(decision_path) as fh:
            obj = json.load(fh)

        rows = obj.get("new_consensus_by_species", {}).get(species, [])
        for row in rows:
            exons = _normalize_exons(row.get("exons", []))
            if not exons:
                continue

            loci.append(
                NewLocusModel(
                    species=species,
                    seqid=row["seqid"],
                    start=min(x[0] for x in exons),
                    end=max(x[1] for x in exons),
                    strand=row["strand"],
                    exons=exons,
                    support_count=_safe_int(row.get("support_count"), 0),
                    total_chain_score=_safe_float(row.get("total_chain_score")),
                    mean_exon_recovery=_safe_float(row.get("mean_exon_recovery")),
                    source_transcripts=list(row.get("source_transcripts", [])),
                )
            )

    return _deduplicate_new_loci(loci)


def write_final_species_gff3s(
    output_dir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_list: list[str],
) -> None:
    """
    For each species, write:
      - final_gff3/<species>.new_loci.gff3
      - final_gff3/<species>.final.gff3
    """
    output_dir_p = Path(output_dir)
    rounds_dir = output_dir_p / "rounds"
    final_dir = output_dir_p / "final_gff3"
    final_dir.mkdir(parents=True, exist_ok=True)

    for species in species_list:
        original_gff3 = Path(annotation_dir) / f"{species}{annotation_suffix}"
        new_loci = load_all_new_loci_for_species(rounds_dir, species)

        write_new_loci_gff3(final_dir / f"{species}.new_loci.gff3", species, new_loci)
        merge_original_and_new_gff3(
            original_gff3=original_gff3,
            final_gff3=final_dir / f"{species}.final.gff3",
            species=species,
            new_loci=new_loci,
        )
