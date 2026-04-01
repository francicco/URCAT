from __future__ import annotations


def interval_overlap_fraction(proj_start: int, proj_end: int, locus_start: int, locus_end: int) -> float:
    """1-based inclusive overlap fraction of the projected span covered by the locus."""
    overlap_start = max(proj_start, locus_start)
    overlap_end = min(proj_end, locus_end)
    overlap = max(0, overlap_end - overlap_start + 1)
    proj_len = max(1, proj_end - proj_start + 1)
    return overlap / proj_len


def strand_consistency(projected, locus) -> float:
    return 1.0 if getattr(projected, "strand", None) == getattr(locus, "strand", None) else 0.0


def orientation_consistency(projected, locus=None) -> float:
    if locus is not None:
        return strand_consistency(projected, locus)
    return 1.0 if getattr(projected, "chain_orientation", None) == "forward" else 0.7
