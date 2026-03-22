from __future__ import annotations

from collections import defaultdict
from statistics import mean

from comparative_annotator.workflow.fragmented_models import (
    FragmentModel,
    LogicalProjectedLocus,
    ProjectedBlock,
)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    xs = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [xs[0]]
    for s, e in xs[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def summarize_projected_blocks(
    target_species: str,
    source_species: str | None,
    source_transcripts: list[str],
    blocks: list[ProjectedBlock],
) -> LogicalProjectedLocus:
    by_seqid: dict[str, list[ProjectedBlock]] = defaultdict(list)
    for b in blocks:
        by_seqid[b.target_seqid].append(b)

    fragments: list[FragmentModel] = []
    seqid_bp: dict[str, int] = {}

    for seqid, seq_blocks in sorted(by_seqid.items()):
        strand_counts: dict[str, int] = defaultdict(int)
        for b in seq_blocks:
            strand_counts[b.target_strand] += 1
        best_strand = sorted(strand_counts.items(), key=lambda x: (-x[1], x[0]))[0][0]

        exons = merge_intervals(
            [(b.target_start, b.target_end) for b in seq_blocks if b.target_strand == best_strand]
        )
        if not exons:
            exons = merge_intervals([(b.target_start, b.target_end) for b in seq_blocks])

        frag = FragmentModel(seqid=seqid, strand=best_strand, exons=exons)
        fragments.append(frag)
        seqid_bp[seqid] = frag.exon_bp

    total_bp = sum(seqid_bp.values())
    dominant_seqid = sorted(seqid_bp.items(), key=lambda x: (-x[1], x[0]))[0][0]
    dominant_bp_fraction = seqid_bp[dominant_seqid] / max(1, total_bp)

    chain_scores = [b.chain_score for b in blocks if b.chain_score is not None]
    exon_recoveries = [b.exon_recovery for b in blocks if b.exon_recovery is not None]

    is_fragmented = len(seqid_bp) > 1
    locus_class, locus_status, notes = classify_projected_locus(
        n_target_seqids=len(seqid_bp),
        dominant_bp_fraction=dominant_bp_fraction,
        mean_exon_recovery=(mean(exon_recoveries) if exon_recoveries else 0.0),
    )

    return LogicalProjectedLocus(
        species=target_species,
        source_species=source_species,
        source_transcripts=sorted(set(source_transcripts)),
        support_count=len(set(source_transcripts)),
        total_chain_score=sum(chain_scores) if chain_scores else 0.0,
        mean_exon_recovery=(mean(exon_recoveries) if exon_recoveries else 0.0),
        fragments=sorted(fragments, key=lambda f: (f.seqid, f.start, f.end)),
        dominant_seqid=dominant_seqid,
        dominant_bp_fraction=dominant_bp_fraction,
        n_target_seqids=len(seqid_bp),
        is_fragmented_across_seqids=is_fragmented,
        locus_class=locus_class,
        locus_status=locus_status,
        notes=notes,
    )


def classify_projected_locus(
    n_target_seqids: int,
    dominant_bp_fraction: float,
    mean_exon_recovery: float,
) -> tuple[str, str, list[str]]:
    notes: list[str] = []

    if n_target_seqids > 1:
        notes.append("projected_exons_split_across_multiple_target_seqids")
        if dominant_bp_fraction < 0.8:
            notes.append("no_single_dominant_seqid")
        return "AF", "fragmented_locus", notes

    if mean_exon_recovery >= 0.9:
        return "I", "new_locus", notes

    if mean_exon_recovery >= 0.4:
        notes.append("partial_exon_recovery")
        return "PI", "partial_locus", notes

    notes.append("weak_projection_support")
    return "UL", "uncertain_locus", notes
