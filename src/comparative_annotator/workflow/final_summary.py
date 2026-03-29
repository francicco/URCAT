from __future__ import annotations

import json
from pathlib import Path


RETAINED_FRAGMENT_CLASSES = {
    "split_fragment_high_confidence",
    "split_fragment_structural",
    "split_fragment_ambiguous",
    "partial_fragment_supported",
}


def _read_json(path: str | Path):
    with open(path) as fh:
        return json.load(fh)


def _iter_noncomment_gff3_lines(path: str | Path):
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            yield line


def _count_gff3_features(path: str | Path) -> dict[str, int]:
    counts = {
        "gene": 0,
        "mRNA": 0,
        "transcript": 0,
    }

    for line in _iter_noncomment_gff3_lines(path):
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        feature = parts[2]
        if feature in counts:
            counts[feature] += 1

    return {
        "genes": counts["gene"],
        "transcripts": counts["mRNA"] + counts["transcript"],
    }


def _consensus_key(species: str, item: dict) -> tuple:
    exons = tuple((int(a), int(b)) for a, b in item.get("exons", []))
    src_txs = tuple(sorted(item.get("source_transcripts", [])))
    return (
        species,
        item.get("seqid"),
        item.get("strand"),
        exons,
        src_txs,
    )


def _fragment_key(target_species: str, item: dict) -> tuple:
    target_seqids = tuple(sorted(str(item.get("target_seqids", "")).split(",")))
    recovered = tuple(
        int(x) for x in str(item.get("recovered_exon_numbers", "")).split(",") if x
    )
    return (
        target_species,
        item.get("source_species"),
        item.get("source_transcript"),
        item.get("source_seqid"),
        item.get("source_strand"),
        int(item.get("n_source_exons", 0) or 0),
        target_seqids,
        recovered,
    )


def _collect_new_consensus_counts(output_dir: str | Path) -> dict[str, int]:
    output_dir = Path(output_dir)
    seen_by_species: dict[str, set] = {}

    for decision_path in sorted(output_dir.glob("rounds/round_*/post_round_decision.json")):
        data = _read_json(decision_path)
        for species, items in (data.get("new_consensus_by_species") or {}).items():
            seen_by_species.setdefault(species, set())
            for item in items:
                seen_by_species[species].add(_consensus_key(species, item))

    return {species: len(keys) for species, keys in seen_by_species.items()}


def _collect_retained_fragment_counts(output_dir: str | Path) -> dict[str, dict[str, int]]:
    output_dir = Path(output_dir)
    seen_by_species: dict[str, dict[str, set]] = {}

    for decision_path in sorted(output_dir.glob("rounds/round_*/post_round_decision.json")):
        data = _read_json(decision_path)
        for species, items in (data.get("fragmented_candidates_by_species") or {}).items():
            block = seen_by_species.setdefault(
                species,
                {
                    "all": set(),
                    "multi_seqid": set(),
                    "partial": set(),
                },
            )

            for item in items:
                fragment_class = item.get("fragment_class")
                if fragment_class not in RETAINED_FRAGMENT_CLASSES:
                    continue

                key = _fragment_key(species, item)
                block["all"].add(key)

                n_target_seqids = int(item.get("n_target_seqids", 0) or 0)
                if n_target_seqids >= 2:
                    block["multi_seqid"].add(key)
                else:
                    block["partial"].add(key)

    out = {}
    for species, block in seen_by_species.items():
        out[species] = {
            "retained_fragmented": len(block["all"]),
            "retained_fragmented_multi_seqid": len(block["multi_seqid"]),
            "retained_fragmented_partial": len(block["partial"]),
        }
    return out


def write_final_summary_tsv(
    output_dir: str,
    annotation_paths: dict[str, str],
    species_list: list[str],
) -> str:
    """
    Write a compact final summary table per species.

    Columns:
      species
      orig_genes
      orig_transcripts
      final_genes
      final_transcripts
      new_complete
      retained_fragmented
      retained_fragmented_multi_seqid
      retained_fragmented_partial
      total_new
      net_gene_gain
    """
    output_dir = Path(output_dir)
    final_gff3_dir = output_dir / "final_gff3"
    final_gff3_dir.mkdir(parents=True, exist_ok=True)

    complete_counts = _collect_new_consensus_counts(output_dir)
    fragmented_counts = _collect_retained_fragment_counts(output_dir)

    rows = []
    for species in species_list:
        orig_path = annotation_paths.get(species, "")
        final_path = final_gff3_dir / f"{species}.final.gff3"

        if orig_path and Path(orig_path).exists():
            orig_counts = _count_gff3_features(orig_path)
        else:
            orig_counts = {"genes": 0, "transcripts": 0}

        if final_path.exists():
            final_counts = _count_gff3_features(final_path)
        else:
            final_counts = {"genes": 0, "transcripts": 0}

        new_complete = complete_counts.get(species, 0)

        frag = fragmented_counts.get(
            species,
            {
                "retained_fragmented": 0,
                "retained_fragmented_multi_seqid": 0,
                "retained_fragmented_partial": 0,
            },
        )

        total_new = new_complete + frag["retained_fragmented"]
        net_gene_gain = final_counts["genes"] - orig_counts["genes"]

        rows.append(
            {
                "species": species,
                "orig_genes": orig_counts["genes"],
                "orig_transcripts": orig_counts["transcripts"],
                "final_genes": final_counts["genes"],
                "final_transcripts": final_counts["transcripts"],
                "new_complete": new_complete,
                "retained_fragmented": frag["retained_fragmented"],
                "retained_fragmented_multi_seqid": frag["retained_fragmented_multi_seqid"],
                "retained_fragmented_partial": frag["retained_fragmented_partial"],
                "total_new": total_new,
                "net_gene_gain": net_gene_gain,
            }
        )

    out_path = final_gff3_dir / "summary_final.tsv"
    header = [
        "species",
        "orig_genes",
        "orig_transcripts",
        "final_genes",
        "final_transcripts",
        "new_complete",
        "retained_fragmented",
        "retained_fragmented_multi_seqid",
        "retained_fragmented_partial",
        "total_new",
        "net_gene_gain",
    ]

    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[col]) for col in header) + "\n")

    return str(out_path)
