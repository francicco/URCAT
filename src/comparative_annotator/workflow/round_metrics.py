from __future__ import annotations

from pathlib import Path

from comparative_annotator.workflow.reporting import read_json, write_tsv


ROUND_METRIC_COLUMNS = [
    "round_id",
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
    "n_new_consensus",
    "n_orphan_loci",
    "n_pending_seeds",
    "next_reference_species",
    "stop",
]


def build_round_metrics_rows(round_dir: str | Path) -> list[dict]:
    round_dir = Path(round_dir)

    summary_path = round_dir / "summary.tsv"
    decision_path = round_dir / "post_round_decision.json"

    if not summary_path.exists():
        return []

    decision = read_json(decision_path) if decision_path.exists() else {}

    new_consensus_counts = {
        sp: len(v) for sp, v in (decision.get("new_consensus_by_species") or {}).items()
    }
    orphan_counts = {
        sp: len(v) for sp, v in (decision.get("orphan_loci_by_species") or {}).items()
    }
    pending_counts = {
        sp: len(v) for sp, v in (decision.get("pending_seeds_by_species") or {}).items()
    }

    rows = []
    with open(summary_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            row = dict(zip(header, fields))
            target_species = row["target_species"]

            edge_path = round_dir / f"ref_{row['reference_species']}" / f"target_{target_species}" / "edge_evidence.json"
            n_edges_total = 0
            n_edges_accepted = 0
            if edge_path.exists():
                x = read_json(edge_path)
                n_edges_total = x.get("n_edges", 0)
                n_edges_accepted = sum(1 for e in x.get("edges", []) if e.get("accepted"))

            rows.append(
                {
                    "round_id": int(row["round_id"]),
                    "reference_species": row["reference_species"],
                    "target_species": target_species,
                    "n_processed": int(row["n_processed"]),
                    "n_ok": int(row["n_ok"]),
                    "n_error": int(row["n_error"]),
                    "n_primary": int(row["n_primary"]),
                    "n_alternative_only": int(row["n_alternative_only"]),
                    "n_missing": int(row["n_missing"]),
                    "n_strand_conflict": int(row["n_strand_conflict"]),
                    "n_edges_total": n_edges_total,
                    "n_edges_accepted": n_edges_accepted,
                    "n_new_consensus": new_consensus_counts.get(target_species, 0),
                    "n_orphan_loci": orphan_counts.get(target_species, 0),
                    "n_pending_seeds": pending_counts.get(target_species, 0),
                    "next_reference_species": row["next_reference_species"] or None,
                    "stop": row["stop"],
                }
            )

    return rows


def write_round_metrics_table(round_dir: str | Path) -> str:
    round_dir = Path(round_dir)
    rows = build_round_metrics_rows(round_dir)
    out_path = round_dir / "round_metrics.tsv"
    write_tsv(out_path, rows, columns=ROUND_METRIC_COLUMNS)
    return str(out_path)
