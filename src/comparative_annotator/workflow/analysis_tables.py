from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from comparative_annotator.workflow.reporting import read_json, write_tsv


ROUND_OVERVIEW_COLUMNS = [
    "round_id",
    "seed_species",
    "reference_species",
    "target_species",
    "n_processed",
    "n_ok",
    "n_error",
    "n_primary",
    "n_alternative_only",
    "n_missing",
    "n_strand_conflict",
    "n_edges_total",
    "n_edges_accepted",
    "n_edge_class_A",
    "n_edge_class_B",
    "n_edge_class_C",
    "n_edge_class_D",
    "n_edge_class_E",
    "n_edge_class_N",
    "n_edge_class_P",
    "n_edge_class_X",
    "n_new_consensus",
    "n_orphan_loci",
    "n_pending_seeds",
    "next_reference_species",
    "stop",
]

PROJECTION_CANDIDATE_COLUMNS = [
    "round_id",
    "reference_species",
    "target_species",
    "source_species",
    "source_transcript_id",
    "source_locus_id",
    "seed_kind",
    "status",
    "primary_target_locus",
    "alternative_target_loci",
    "missing_intervals",
    "strand_conflicts",
    "n_primary_targets",
    "n_alternative_targets",
    "n_missing_intervals",
]

EDGE_SUMMARY_COLUMNS = [
    "round_id",
    "reference_species",
    "target_species",
    "source_species",
    "source_locus_id",
    "target_locus_id",
    "edge_origin",
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
    "proj_projected_seqid_match",
    "syn_synteny_class",
    "syn_anchor_ok",
    "syn_anchor_distance",
    "syn_left_neighbor_match",
    "syn_right_neighbor_match",
    "syn_neighbor_hits",
    "syn_neighbor_jaccard",
    "syn_order_score",
    "syn_orientation_score",
    "arch_shared_introns",
    "arch_intron_recall",
    "arch_intron_precision",
    "arch_exon_count_compatible",
    "arch_utr_support",
    "flag_strand_conflict",
    "flag_has_competitor",
]

NOVEL_AND_PENDING_COLUMNS = [
    "round_id",
    "reference_species",
    "species",
    "event_type",
    "selected_for_next_round",
    "seqid",
    "start",
    "end",
    "strand",
    "support_count",
    "total_chain_score",
    "mean_exon_recovery",
    "source_species",
    "source_transcripts",
    "source_loci",
    "best_source_transcript",
    "best_target_transcript",
    "protein_identity",
    "protein_coverage",
    "bitscore",
    "supporting_target_locus",
    "supporting_edge_class",
    "supporting_edge_confidence",
    "locus_id",
    "transcripts",
]


def _infer_source_locus_id(source_transcript_id: str | None) -> str | None:
    if source_transcript_id is None:
        return None
    tx = source_transcript_id.removeprefix("transcript:")
    if "." in tx:
        return tx.rsplit(".", 1)[0]
    return tx


def _safe_int(x: Any, default: int = 0) -> int:
    if x in (None, "", "None"):
        return default
    return int(x)


def _safe_boolish(x: Any) -> Any:
    if x in ("True", True):
        return True
    if x in ("False", False):
        return False
    return x


def _best_edge_score(edge: dict[str, Any]) -> tuple[float, float, float]:
    bitscore = edge.get("seq_bitscore")
    prot_id = edge.get("seq_prot_id")
    prot_cov = edge.get("seq_prot_cov")
    return (
        -1.0 if bitscore is None else float(bitscore),
        -1.0 if prot_id is None else float(prot_id),
        -1.0 if prot_cov is None else float(prot_cov),
    )


def _choose_best_edge(edges: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not edges:
        return None
    return max(edges, key=_best_edge_score)


def _interval_fields(item: dict) -> dict[str, Any]:
    return {
        "seqid": item.get("seqid"),
        "start": item.get("start"),
        "end": item.get("end"),
        "strand": item.get("strand"),
    }


def _load_summary_rows(round_dir: str | Path) -> list[dict[str, Any]]:
    round_dir = Path(round_dir)
    summary_path = round_dir / "summary.tsv"
    if not summary_path.exists():
        return []

    rows: list[dict[str, Any]] = []
    with open(summary_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, vals)))
    return rows


def _build_edge_rows_for_target(target_dir: str | Path) -> list[dict[str, Any]]:
    target_dir = Path(target_dir)
    merged_path = target_dir / "merged.json"
    edge_path = target_dir / "edge_evidence.json"

    if not merged_path.exists() or not edge_path.exists():
        return []

    merged = read_json(merged_path)
    edge_payload = read_json(edge_path)
    rows: list[dict[str, Any]] = []

    for e in edge_payload.get("edges", []):
        rows.append(
            {
                "round_id": merged.get("round_id"),
                "reference_species": merged.get("reference_species"),
                "target_species": merged.get("target_species"),
                "source_species": e.get("source_species"),
                "source_locus_id": e.get("source_locus_id"),
                "target_locus_id": e.get("target_locus_id"),
                "edge_origin": e.get("edge_origin"),
                "accepted": e.get("accepted"),
                "edge_class": e.get("edge_class"),
                "edge_confidence": e.get("edge_confidence"),
                "seq_best_source_tx": e.get("seq_best_source_tx"),
                "seq_best_target_tx": e.get("seq_best_target_tx"),
                "seq_prot_id": e.get("seq_prot_id"),
                "seq_prot_cov": e.get("seq_prot_cov"),
                "seq_bitscore": e.get("seq_bitscore"),
                "proj_cov": e.get("proj_cov"),
                "proj_same_strand": e.get("proj_same_strand"),
                "proj_anchor_distance": e.get("proj_anchor_distance"),
                "proj_projected_seqid_match": e.get("proj_projected_seqid_match"),
                "syn_synteny_class": e.get("syn_synteny_class"),
                "syn_anchor_ok": e.get("syn_anchor_ok"),
                "syn_anchor_distance": e.get("syn_anchor_distance"),
                "syn_left_neighbor_match": e.get("syn_left_neighbor_match"),
                "syn_right_neighbor_match": e.get("syn_right_neighbor_match"),
                "syn_neighbor_hits": e.get("syn_neighbor_hits"),
                "syn_neighbor_jaccard": e.get("syn_neighbor_jaccard"),
                "syn_order_score": e.get("syn_order_score"),
                "syn_orientation_score": e.get("syn_orientation_score"),
                "arch_shared_introns": e.get("arch_shared_introns"),
                "arch_intron_recall": e.get("arch_intron_recall"),
                "arch_intron_precision": e.get("arch_intron_precision"),
                "arch_exon_count_compatible": e.get("arch_exon_count_compatible"),
                "arch_utr_support": e.get("arch_utr_support"),
                "flag_strand_conflict": e.get("flag_strand_conflict"),
                "flag_has_competitor": e.get("flag_has_competitor"),
            }
        )

    return rows


def _build_projection_candidate_rows_for_target(target_dir: str | Path) -> list[dict[str, Any]]:
    target_dir = Path(target_dir)
    merged_path = target_dir / "merged.json"

    if not merged_path.exists():
        return []

    merged = read_json(merged_path)
    target_species = merged.get("target_species")
    rows: list[dict[str, Any]] = []

    for r in merged.get("results", []):
        if r.get("status") != "ok":
            continue

        primary_map = r.get("primary") or {}
        alt_map = r.get("alternatives") or {}
        missing_map = r.get("missing_annotations") or {}
        strand_map = r.get("strand_conflicts") or {}

        primary_target = primary_map.get(target_species)
        alt_targets = alt_map.get(target_species, [])
        missing_intervals = missing_map.get(target_species, [])
        strand_conflicts = strand_map.get(target_species, [])

        if primary_target is not None or alt_targets:
            status = "projected"
        elif missing_intervals:
            status = "missing"
        else:
            status = "no_hit"

        source_tx = r.get("source_transcript")
        rows.append(
            {
                "round_id": merged.get("round_id"),
                "reference_species": merged.get("reference_species"),
                "target_species": target_species,
                "source_species": r.get("source_species"),
                "source_transcript_id": source_tx,
                "source_locus_id": _infer_source_locus_id(source_tx),
                "seed_kind": r.get("seed_kind"),
                "status": status,
                "primary_target_locus": primary_target,
                "alternative_target_loci": alt_targets,
                "missing_intervals": missing_intervals,
                "strand_conflicts": strand_conflicts,
                "n_primary_targets": 0 if primary_target is None else 1,
                "n_alternative_targets": len(alt_targets),
                "n_missing_intervals": len(missing_intervals),
            }
        )

    return rows


def _build_best_edge_support_index(round_dir: str | Path) -> dict[tuple[str, str], dict[str, Any]]:
    round_dir = Path(round_dir)
    best_by_species_and_tx: dict[tuple[str, str], list[dict[str, Any]]] = {}

    for edge_json in sorted(round_dir.rglob("edge_evidence.json")):
        payload = read_json(edge_json)
        target_species = payload.get("target_species")

        for edge in payload.get("edges", []):
            source_tx = edge.get("seq_best_source_tx")
            if not source_tx:
                continue
            key = (target_species, source_tx)
            best_by_species_and_tx.setdefault(key, []).append(edge)

    out: dict[tuple[str, str], dict[str, Any]] = {}
    for key, edges in best_by_species_and_tx.items():
        best = _choose_best_edge(edges)
        if best is not None:
            out[key] = best
    return out


def _best_support_for_row(
    row_species: str,
    source_transcripts: list[str] | None,
    support_index: dict[tuple[str, str], dict[str, Any]],
) -> dict[str, Any]:
    if not source_transcripts:
        return {
            "source_species": None,
            "best_source_transcript": None,
            "best_target_transcript": None,
            "protein_identity": None,
            "protein_coverage": None,
            "bitscore": None,
            "supporting_target_locus": None,
            "supporting_edge_class": None,
            "supporting_edge_confidence": None,
        }

    candidate_edges = []
    for tx in source_transcripts:
        edge = support_index.get((row_species, tx))
        if edge is not None:
            candidate_edges.append(edge)

    best = _choose_best_edge(candidate_edges)
    if best is None:
        return {
            "source_species": None,
            "best_source_transcript": None,
            "best_target_transcript": None,
            "protein_identity": None,
            "protein_coverage": None,
            "bitscore": None,
            "supporting_target_locus": None,
            "supporting_edge_class": None,
            "supporting_edge_confidence": None,
        }

    return {
        "source_species": best.get("source_species"),
        "best_source_transcript": best.get("seq_best_source_tx"),
        "best_target_transcript": best.get("seq_best_target_tx"),
        "protein_identity": best.get("seq_prot_id"),
        "protein_coverage": best.get("seq_prot_cov"),
        "bitscore": best.get("seq_bitscore"),
        "supporting_target_locus": best.get("target_locus_id"),
        "supporting_edge_class": best.get("edge_class"),
        "supporting_edge_confidence": best.get("edge_confidence"),
    }


def _build_novel_and_pending_rows(round_dir: str | Path) -> list[dict[str, Any]]:
    round_dir = Path(round_dir)
    decision_path = round_dir / "post_round_decision.json"
    if not decision_path.exists():
        return []

    support_index = _build_best_edge_support_index(round_dir)
    x = read_json(decision_path)
    round_id = x.get("round_id")
    reference_species = x.get("reference_species")
    next_reference_species = x.get("next_reference_species")

    rows: list[dict[str, Any]] = []

    for species, items in (x.get("new_consensus_by_species") or {}).items():
        for item in items:
            support = _best_support_for_row(species, item.get("source_transcripts"), support_index)
            row = {
                "round_id": round_id,
                "reference_species": reference_species,
                "species": species,
                "event_type": "new_consensus",
                "selected_for_next_round": species == next_reference_species,
                "support_count": item.get("support_count"),
                "total_chain_score": item.get("total_chain_score"),
                "mean_exon_recovery": item.get("mean_exon_recovery"),
                "source_transcripts": item.get("source_transcripts"),
                "source_loci": item.get("source_loci"),
                "locus_id": item.get("locus_id"),
                "transcripts": item.get("transcripts"),
            }
            row.update(_interval_fields(item))
            row.update(support)
            rows.append(row)

    for species, items in (x.get("orphan_loci_by_species") or {}).items():
        for item in items:
            row = {
                "round_id": round_id,
                "reference_species": reference_species,
                "species": species,
                "event_type": "orphan_native_locus",
                "selected_for_next_round": species == next_reference_species,
                "support_count": None,
                "total_chain_score": None,
                "mean_exon_recovery": None,
                "source_transcripts": None,
                "source_loci": None,
                "source_species": None,
                "best_source_transcript": None,
                "best_target_transcript": None,
                "protein_identity": None,
                "protein_coverage": None,
                "bitscore": None,
                "supporting_target_locus": None,
                "supporting_edge_class": None,
                "supporting_edge_confidence": None,
                "locus_id": item.get("locus_id"),
                "transcripts": item.get("transcripts"),
            }
            row.update(_interval_fields(item))
            rows.append(row)

    pending_map = x.get("pending_seeds_by_species") or x.get("pending_frontiers_by_species") or {}
    for species, items in pending_map.items():
        for item in items:
            support = _best_support_for_row(species, item.get("source_transcripts"), support_index)
            row = {
                "round_id": round_id,
                "reference_species": reference_species,
                "species": species,
                "event_type": "pending_seed",
                "selected_for_next_round": species == next_reference_species,
                "support_count": item.get("support_count"),
                "total_chain_score": item.get("total_chain_score"),
                "mean_exon_recovery": item.get("mean_exon_recovery"),
                "source_transcripts": item.get("source_transcripts"),
                "source_loci": item.get("source_loci"),
                "locus_id": item.get("locus_id"),
                "transcripts": item.get("transcripts"),
            }
            row.update(_interval_fields(item))
            row.update(support)
            rows.append(row)

    return rows


def _build_round_overview_rows(round_dir: str | Path) -> list[dict[str, Any]]:
    round_dir = Path(round_dir)
    summary_rows = _load_summary_rows(round_dir)

    decision_path = round_dir / "post_round_decision.json"
    decision = read_json(decision_path) if decision_path.exists() else {}

    new_consensus_counts = {
        sp: len(v) for sp, v in (decision.get("new_consensus_by_species") or {}).items()
    }
    orphan_counts = {
        sp: len(v) for sp, v in (decision.get("orphan_loci_by_species") or {}).items()
    }
    pending_counts = {
        sp: len(v)
        for sp, v in (
            decision.get("pending_seeds_by_species")
            or decision.get("pending_frontiers_by_species")
            or {}
        ).items()
    }

    out_rows = []

    for srow in summary_rows:
        reference_species = srow["reference_species"]
        target_species = srow["target_species"]

        edge_path = round_dir / f"ref_{reference_species}" / f"target_{target_species}" / "edge_evidence.json"
        class_counts: Counter[str] = Counter()
        n_edges_total = 0
        n_edges_accepted = 0

        if edge_path.exists():
            edge_payload = read_json(edge_path)
            edge_rows = edge_payload.get("edges", [])
            n_edges_total = edge_payload.get("n_edges", len(edge_rows))
            n_edges_accepted = sum(1 for e in edge_rows if e.get("accepted") is True)
            class_counts = Counter(e.get("edge_class") for e in edge_rows if e.get("edge_class"))

        out_rows.append({
            "round_id": srow.get("round_id"),
            "seed_species": srow.get("seed_species"),
            "reference_species": reference_species,
            "target_species": target_species,
            "n_processed": srow.get("n_processed"),
            "n_ok": srow.get("n_ok"),
            "n_error": srow.get("n_error"),
            "n_primary": srow.get("n_primary"),
            "n_alternative_only": srow.get("n_alternative_only"),
            "n_missing": srow.get("n_missing"),
            "n_strand_conflict": srow.get("n_strand_conflict"),
            "n_edges_total": n_edges_total,
            "n_edges_accepted": n_edges_accepted,
            "n_edge_class_A": class_counts.get("A", 0),
            "n_edge_class_B": class_counts.get("B", 0),
            "n_edge_class_C": class_counts.get("C", 0),
            "n_edge_class_D": class_counts.get("D", 0),
            "n_edge_class_E": class_counts.get("E", 0),
            "n_edge_class_N": class_counts.get("N", 0),
            "n_edge_class_P": class_counts.get("P", 0),
            "n_edge_class_X": class_counts.get("X", 0),
            "n_new_consensus": new_consensus_counts.get(target_species, 0),
            "n_orphan_loci": orphan_counts.get(target_species, 0),
            "n_pending_seeds": pending_counts.get(target_species, 0),
            "next_reference_species": srow.get("next_reference_species"),
            "stop": srow.get("stop"),
        })

    return out_rows


def write_edge_summary_tables_for_target(target_dir: str | Path) -> tuple[str, str]:
    target_dir = Path(target_dir)

    edge_rows = _build_edge_rows_for_target(target_dir)
    edge_summary_path = target_dir / "edge_summary.tsv"
    accepted_edges_path = target_dir / "accepted_edges.tsv"

    write_tsv(edge_summary_path, edge_rows, columns=EDGE_SUMMARY_COLUMNS)
    write_tsv(
        accepted_edges_path,
        [r for r in edge_rows if r.get("accepted") is True],
        columns=EDGE_SUMMARY_COLUMNS,
    )
    return str(edge_summary_path), str(accepted_edges_path)


def write_projection_candidates_table_for_target(target_dir: str | Path) -> str:
    target_dir = Path(target_dir)
    rows = _build_projection_candidate_rows_for_target(target_dir)
    out_path = target_dir / "projection_candidates.tsv"
    write_tsv(out_path, rows, columns=PROJECTION_CANDIDATE_COLUMNS)
    return str(out_path)


def write_round_overview_table(round_dir: str | Path) -> str:
    round_dir = Path(round_dir)
    rows = _build_round_overview_rows(round_dir)
    out_path = round_dir / "round_overview.tsv"
    write_tsv(out_path, rows, columns=ROUND_OVERVIEW_COLUMNS)
    return str(out_path)


def write_novel_and_pending_table(round_dir: str | Path) -> str:
    round_dir = Path(round_dir)
    rows = _build_novel_and_pending_rows(round_dir)
    out_path = round_dir / "novel_and_pending.tsv"
    write_tsv(out_path, rows, columns=NOVEL_AND_PENDING_COLUMNS)
    return str(out_path)


def write_all_analysis_tables_for_target(target_dir: str | Path) -> None:
    write_projection_candidates_table_for_target(target_dir)
    write_edge_summary_tables_for_target(target_dir)


def write_all_analysis_tables_for_round(round_dir: str | Path) -> None:
    round_dir = Path(round_dir)
    write_round_overview_table(round_dir)
    write_novel_and_pending_table(round_dir)

    for target_dir in sorted(round_dir.rglob("target_*")):
        if target_dir.is_dir():
            write_all_analysis_tables_for_target(target_dir)
