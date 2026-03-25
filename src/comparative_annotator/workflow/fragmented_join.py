from __future__ import annotations

from collections import defaultdict

from comparative_annotator.models.fragmented_locus import (
    FragmentedBlock,
    FragmentedComparativeLocus,
)


def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else 0.0


def build_fragmented_candidate(
    *,
    source_species: str,
    source_transcript: str,
    target_species: str,
    source_seqid: str,
    source_strand: str,
    n_source_exons: int,
    projected_blocks: list[dict],
    min_exon_recovery_fraction: float = 0.35,
    min_recovered_exons: int = 2,
    max_target_seqids: int = 4,
) -> FragmentedComparativeLocus | None:
    """
    projected_blocks rows must contain:
      source_exon_number, target_seqid, target_start, target_end, target_strand, chain_score
    """

    if not projected_blocks:
        return None

    by_exon = defaultdict(list)
    for row in projected_blocks:
        by_exon[row["source_exon_number"]].append(row)

    chosen = []
    for exon_no in sorted(by_exon):
        rows = by_exon[exon_no]
        rows = sorted(rows, key=lambda x: (x.get("chain_score") or 0.0), reverse=True)
        chosen.append(rows[0])

    recovered_exon_numbers = [r["source_exon_number"] for r in chosen]
    n_recovered = len(recovered_exon_numbers)
    exon_recovery_fraction = n_recovered / max(1, n_source_exons)

    target_seqids = []
    for r in chosen:
        if r["target_seqid"] not in target_seqids:
            target_seqids.append(r["target_seqid"])

    if n_recovered < min_recovered_exons:
        return None
    if exon_recovery_fraction < min_exon_recovery_fraction:
        return None
    if len(target_seqids) < 2:
        return None
    if len(target_seqids) > max_target_seqids:
        return None

    strand_set = {r["target_strand"] for r in chosen}
    strand_consistent = len(strand_set) == 1

    exon_order_consistent = True
    per_seqid = defaultdict(list)
    for r in chosen:
        per_seqid[r["target_seqid"]].append(r)

    for seqid, rows in per_seqid.items():
        rows = sorted(rows, key=lambda x: x["source_exon_number"])
        coords = [r["target_start"] for r in rows]
        if coords != sorted(coords):
            exon_order_consistent = False
            break

    chain_scores = [(r.get("chain_score") or 0.0) for r in chosen]
    total_chain_score = sum(chain_scores)

    # Simple scoring function; tune later.
    join_score = (
        4.0 * exon_recovery_fraction
        + 1.5 * (1.0 if strand_consistent else 0.0)
        + 1.5 * (1.0 if exon_order_consistent else 0.0)
        + min(2.0, _mean(chain_scores) / 1000.0)
        - 0.75 * (len(target_seqids) - 1)
    )

    if join_score < 2.5:
        return None

    return FragmentedComparativeLocus(
        source_species=source_species,
        source_transcript=source_transcript,
        target_species=target_species,
        source_seqid=source_seqid,
        source_strand=source_strand,
        n_source_exons=n_source_exons,
        target_seqids=target_seqids,
        blocks=[
            FragmentedBlock(
                source_exon_number=r["source_exon_number"],
                target_seqid=r["target_seqid"],
                target_start=r["target_start"],
                target_end=r["target_end"],
                target_strand=r["target_strand"],
                chain_score=r.get("chain_score"),
            )
            for r in chosen
        ],
        recovered_exon_numbers=recovered_exon_numbers,
        exon_recovery_fraction=exon_recovery_fraction,
        strand_consistent=strand_consistent,
        exon_order_consistent=exon_order_consistent,
        total_chain_score=total_chain_score,
        join_score=join_score,
        status="fragmented",
    )
