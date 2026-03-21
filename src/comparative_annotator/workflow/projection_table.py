from __future__ import annotations

from pathlib import Path

from comparative_annotator.workflow.reporting import read_json, write_tsv


PROJECTION_COLUMNS = [
    "round_id",
    "reference_species",
    "target_species",
    "source_species",
    "source_locus_id",
    "source_transcript_id",
    "seed_kind",
    "target_locus_id",
    "edge_origin",
    "projection_status",
    "primary_target_locus",
    "alternative_target_loci",
    "missing_interval",
    "strand_conflict",
    "accepted",
    "edge_class",
    "edge_confidence",
    "seq_best_source_tx",
    "seq_best_target_tx",
    "seq_prot_id",
    "seq_prot_cov",
    "seq_bitscore",
    "proj_cov",
    "proj_same_strand",
    "proj_anchor_distance",
    "syn_synteny_class",
    "syn_anchor_ok",
    "arch_shared_introns",
    "arch_intron_recall",
    "arch_intron_precision",
    "arch_exon_count_compatible",
]


def _edge_index(edge_rows: list[dict]) -> dict[tuple[str, str, str, str], dict]:
    idx = {}
    for row in edge_rows:
        key = (
            row.get("source_species"),
            row.get("source_locus_id"),
            row.get("target_species"),
            row.get("target_locus_id"),
        )
        idx[key] = row
    return idx


def _infer_source_locus_id(source_transcript_id: str | None) -> str | None:
    if source_transcript_id is None:
        return None
    tx = source_transcript_id.removeprefix("transcript:")
    if "." in tx:
        return tx.rsplit(".", 1)[0]
    return tx


def build_projection_rows_for_target(target_dir: str | Path) -> list[dict]:
    target_dir = Path(target_dir)

    merged_path = target_dir / "merged.json"
    edge_path = target_dir / "edge_evidence.json"

    if not merged_path.exists():
        return []

    merged = read_json(merged_path)
    edges = read_json(edge_path) if edge_path.exists() else {"edges": []}
    edge_idx = _edge_index(edges.get("edges", []))

    rows = []

    for r in merged.get("results", []):
        if r.get("status") != "ok":
            continue

        round_id = merged.get("round_id")
        reference_species = merged.get("reference_species")
        target_species = merged.get("target_species")
        source_species = r.get("source_species")
        source_transcript_id = r.get("source_transcript")
        source_locus_id = _infer_source_locus_id(source_transcript_id)

        primary_target_locus = (r.get("primary") or {}).get(target_species)
        alt_target_loci = (r.get("alternatives") or {}).get(target_species, [])
        missing_tokens = (r.get("missing_annotations") or {}).get(target_species, [])
        strand_conflicts = (r.get("strand_conflicts") or {}).get(target_species, [])

        candidate_target_loci = []
        if primary_target_locus is not None:
            candidate_target_loci.append((primary_target_locus, "primary"))
        for alt in alt_target_loci:
            candidate_target_loci.append((alt, "alternative"))

        if not candidate_target_loci:
            rows.append(
                {
                    "round_id": round_id,
                    "reference_species": reference_species,
                    "target_species": target_species,
                    "source_species": source_species,
                    "source_locus_id": source_locus_id,
                    "source_transcript_id": source_transcript_id,
                    "seed_kind": r.get("seed_kind"),
                    "target_locus_id": None,
                    "edge_origin": None,
                    "projection_status": "missing" if missing_tokens else "no_hit",
                    "primary_target_locus": primary_target_locus,
                    "alternative_target_loci": alt_target_loci,
                    "missing_interval": missing_tokens,
                    "strand_conflict": strand_conflicts,
                    "accepted": None,
                    "edge_class": None,
                    "edge_confidence": None,
                    "seq_best_source_tx": None,
                    "seq_best_target_tx": None,
                    "seq_prot_id": None,
                    "seq_prot_cov": None,
                    "seq_bitscore": None,
                    "proj_cov": None,
                    "proj_same_strand": None,
                    "proj_anchor_distance": None,
                    "syn_synteny_class": None,
                    "syn_anchor_ok": None,
                    "arch_shared_introns": None,
                    "arch_intron_recall": None,
                    "arch_intron_precision": None,
                    "arch_exon_count_compatible": None,
                }
            )
            continue

        for target_locus_id, edge_origin in candidate_target_loci:
            edge = edge_idx.get(
                (source_species, source_locus_id, target_species, target_locus_id),
                {},
            )

            rows.append(
                {
                    "round_id": round_id,
                    "reference_species": reference_species,
                    "target_species": target_species,
                    "source_species": source_species,
                    "source_locus_id": source_locus_id,
                    "source_transcript_id": source_transcript_id,
                    "seed_kind": r.get("seed_kind"),
                    "target_locus_id": target_locus_id,
                    "edge_origin": edge_origin,
                    "projection_status": "projected",
                    "primary_target_locus": primary_target_locus,
                    "alternative_target_loci": alt_target_loci,
                    "missing_interval": missing_tokens,
                    "strand_conflict": strand_conflicts,
                    "accepted": edge.get("accepted"),
                    "edge_class": edge.get("edge_class"),
                    "edge_confidence": edge.get("edge_confidence"),
                    "seq_best_source_tx": edge.get("seq_best_source_tx"),
                    "seq_best_target_tx": edge.get("seq_best_target_tx"),
                    "seq_prot_id": edge.get("seq_prot_id"),
                    "seq_prot_cov": edge.get("seq_prot_cov"),
                    "seq_bitscore": edge.get("seq_bitscore"),
                    "proj_cov": edge.get("proj_cov"),
                    "proj_same_strand": edge.get("proj_same_strand"),
                    "proj_anchor_distance": edge.get("proj_anchor_distance"),
                    "syn_synteny_class": edge.get("syn_synteny_class"),
                    "syn_anchor_ok": edge.get("syn_anchor_ok"),
                    "arch_shared_introns": edge.get("arch_shared_introns"),
                    "arch_intron_recall": edge.get("arch_intron_recall"),
                    "arch_intron_precision": edge.get("arch_intron_precision"),
                    "arch_exon_count_compatible": edge.get("arch_exon_count_compatible"),
                }
            )

    return rows


def write_projection_evidence_table(target_dir: str | Path) -> str:
    target_dir = Path(target_dir)
    rows = build_projection_rows_for_target(target_dir)
    out_path = target_dir / "projection_evidence.tsv"
    write_tsv(out_path, rows, columns=PROJECTION_COLUMNS)
    return str(out_path)
