"""QC reporting utilities for URCAT."""
from __future__ import annotations

def read_json(path: str | Path):
    path = Path(path)
    with open(path) as fh:
        return json.load(fh)


def write_json(path: str | Path, obj) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)

def format_projection_report(
    seed_transcript_id: str,
    seed_species: str,
    target_results: list[dict],
) -> str:
    """
    Format a human-readable projection diagnostic report.

    Parameters
    ----------
    seed_transcript_id:
        e.g. 'transcript:Hmel202001oG3.1'
    seed_species:
        e.g. 'Hmel'
    target_results:
        List of per-target dicts with keys:
          target_species, seqid, start, end, exon_count, orientation,
          coverage, primary, strand_conflicts, missing_annotation,
          nearest_same_strand_locus, nearest_locus_distance
    """
    lines = [f"Seed: {seed_transcript_id} ({seed_species})", ""]
    for r in target_results:
        lines.append(f"Target: {r['target_species']}")
        lines.append("Projected candidate:")
        lines.append(f"  seqid: {r.get('seqid', 'n/a')}")
        lines.append(f"  span: {r.get('start')}-{r.get('end')}")
        lines.append(f"  exon_count: {r.get('exon_count', 0)}")
        lines.append(f"  orientation: {r.get('orientation', 'n/a')}")
        lines.append(f"  coverage: {r.get('coverage', 0)}")
        lines.append("Result:")
        lines.append(f"  primary: {r.get('primary', 'none')}")
        lines.append(f"  strand_conflicts: {r.get('strand_conflicts', 'none')}")
        lines.append(f"  missing_annotation: {r.get('missing_annotation', 'none')}")
        dist = r.get("nearest_locus_distance")
        locus = r.get("nearest_same_strand_locus")
        if locus:
            lines.append(f"  nearest_same_strand_locus: {locus} (distance {dist})")
        lines.append("")
    return "\n".join(lines)
