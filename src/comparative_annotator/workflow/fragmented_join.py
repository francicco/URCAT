from __future__ import annotations

import math
from collections import defaultdict

from comparative_annotator.models.fragmented_locus import (
    FragmentedBlock,
    FragmentedComparativeLocus,
)


def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else 0.0

def _normalize_tx_id(x: str) -> str:
    x = x.split()[0]
    if x.startswith("transcript:"):
        x = x[len("transcript:"):]
    return x


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []

    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


def _intervals_total_length(intervals: list[tuple[int, int]]) -> int:
    return sum((end - start + 1) for start, end in intervals)
    
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
    if not projected_blocks:
        print("FRAGREJECT no_projected_blocks", source_species, source_transcript, target_species)
        return None

    by_exon = defaultdict(list)
    for row in projected_blocks:
        by_exon[row["source_exon_number"]].append(row)

    chosen = []
    for exon_no in sorted(by_exon):
        rows = sorted(
            by_exon[exon_no],
            key=lambda x: (x.get("chain_score") or 0.0),
            reverse=True,
        )
        chosen.append(rows[0])

    recovered_exon_numbers = [r["source_exon_number"] for r in chosen]
    n_recovered = len(recovered_exon_numbers)
    exon_recovery_fraction = n_recovered / max(1, n_source_exons)

    target_seqids = []
    for r in chosen:
        if r["target_seqid"] not in target_seqids:
            target_seqids.append(r["target_seqid"])

    print(
        "FRAGDEBUG",
        source_species,
        source_transcript,
        target_species,
        "n_source_exons=", n_source_exons,
        "n_recovered=", n_recovered,
        "recovery=", round(exon_recovery_fraction, 3),
        "n_target_seqids=", len(target_seqids),
        "target_seqids=", ",".join(target_seqids),
    )

    if n_recovered < min_recovered_exons:
        print("FRAGREJECT too_few_recovered_exons", source_transcript, n_recovered)
        return None

    if exon_recovery_fraction < min_exon_recovery_fraction:
        print("FRAGREJECT low_recovery_fraction", source_transcript, exon_recovery_fraction)
        return None

    multi_seqid_fragment = len(target_seqids) >= 2
    single_seqid_partial = (
        len(target_seqids) == 1
        and n_recovered >= min_recovered_exons
        and exon_recovery_fraction >= min_exon_recovery_fraction
        and exon_recovery_fraction < 0.8
    )

    if not (multi_seqid_fragment or single_seqid_partial):
        print("FRAGREJECT neither_multi_seqid_nor_single_partial", source_transcript)
        return None

    if len(target_seqids) > max_target_seqids:
        print("FRAGREJECT too_many_target_seqids", source_transcript, len(target_seqids))
        return None

    strand_set = {r["target_strand"] for r in chosen}
    strand_consistent = len(strand_set) == 1

    exon_order_consistent = True
    per_seqid = defaultdict(list)
    for r in chosen:
        per_seqid[r["target_seqid"]].append(r)

    for seqid, rows in per_seqid.items():
        rows = sorted(rows, key=lambda x: x["source_exon_number"])
        starts = [r["target_start"] for r in rows]
        if starts != sorted(starts):
            exon_order_consistent = False
            break

    chain_scores = [(r.get("chain_score") or 0.0) for r in chosen]
    total_chain_score = sum(chain_scores)

    join_score = (
        4.0 * exon_recovery_fraction
        + 1.5 * (1.0 if strand_consistent else 0.0)
        + 1.5 * (1.0 if exon_order_consistent else 0.0)
        + min(2.0, _mean(chain_scores) / 1000.0)
        - 0.75 * (len(target_seqids) - 1)
    )

    print(
        "FRAGDEBUG_SCORE",
        source_transcript,
        "strand_consistent=", strand_consistent,
        "exon_order_consistent=", exon_order_consistent,
        "join_score=", round(join_score, 3),
    )

    if join_score < 2.5:
        print("FRAGREJECT low_join_score", source_transcript, join_score)
        return None

    print("FRAGACCEPT", source_transcript, target_species)

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


def add_protein_support_to_candidate(
    candidate,
    diamond_hits_for_source: dict[tuple[str, str], dict],
    source_protein_length: int,
):
    """
    Attach protein support to a fragmented/disrupted candidate.

    Coverage is computed from the union of covered query intervals, not from
    the sum of HSP lengths, to avoid coverage > 1.0 due to overlapping hits.
    """
    candidate.aa_identity_mean = None
    candidate.aa_coverage_fraction = None
    candidate.aa_bitscore = None
    candidate.protein_support_class = "no_protein_evidence"

    if not diamond_hits_for_source:
        return candidate

    source_tx_norm = _normalize_tx_id(candidate.source_transcript)
    target_seqids = set(candidate.target_seqids)

    matched_hits = []

    for (qseqid, sseqid), hit in diamond_hits_for_source.items():
        if _normalize_tx_id(qseqid) != source_tx_norm:
            continue

        # If your sseqid encodes target transcript IDs rather than seqids,
        # you may want to relax this filter later. For now we keep all query-matching hits.
        matched_hits.append(((qseqid, sseqid), hit))

    if not matched_hits:
        return candidate

    # Identity and bitscore: use best-supported hit set
    pidents = []
    bitscores = []

    # Coverage: union of query intervals if present, otherwise fall back safely
    query_intervals = []

    for (_, _), hit in matched_hits:
        pid = hit.get("pid")
        bitscore = hit.get("bitscore")

        if pid is not None:
            pidents.append(float(pid))
        if bitscore is not None:
            bitscores.append(float(bitscore))

        # Preferred: explicit query coordinates from DIAMOND parser
        qstart = hit.get("qstart")
        qend = hit.get("qend")
        if qstart is not None and qend is not None:
            query_intervals.append((int(qstart), int(qend))) / 100.0

    if pidents:
        candidate.aa_identity_mean = sum(pidents) / len(pidents)

    if bitscores:
        candidate.aa_bitscore = max(bitscores)

    if query_intervals:
        merged = _merge_intervals(query_intervals)
        covered_query_len = _intervals_total_length(merged)
        candidate.aa_coverage_fraction = min(covered_query_len / max(1, source_protein_length), 1.0)
    else:
        # Fallback if qstart/qend are not available yet in parsed DIAMOND results.
        # Use best single-hit alignment length only, not sum across hits.
        aln_lens = [
            int(hit["aln_len"])
            for (_, _), hit in matched_hits
            if hit.get("aln_len") is not None
        ]
        if aln_lens:
            candidate.aa_coverage_fraction = min(
                max(aln_lens) / max(1, source_protein_length),
                1.0,
            )

    pid = candidate.aa_identity_mean
    cov = candidate.aa_coverage_fraction

    if pid is None or cov is None:
        candidate.protein_support_class = "no_protein_evidence"
    elif pid >= 0.80 and cov >= 0.80:
        candidate.protein_support_class = "strong"
    elif pid >= 0.60 and cov >= 0.50:
        candidate.protein_support_class = "moderate"
    elif pid >= 0.45 and cov >= 0.30:
        candidate.protein_support_class = "weak"
    else:
        candidate.protein_support_class = "poor"

    return candidate


def classify_fragmented_candidate(
    candidate: FragmentedComparativeLocus,
) -> FragmentedComparativeLocus:
    """
    Final classification using structure + protein support.
    """
    pid = candidate.aa_identity_mean or 0.0
    cov = candidate.aa_coverage_fraction if candidate.aa_coverage_fraction is not None else 0.0

    strong_structure = (
        candidate.exon_recovery_fraction >= 0.40
        and candidate.strand_consistent
        and candidate.exon_order_consistent
        and len(candidate.target_seqids) <= 3
        and candidate.join_score >= 3.0
    )

    moderate_structure = (
        candidate.exon_recovery_fraction >= 0.25
        and candidate.join_score >= 2.5
    )

    strong_protein = pid >= 0.80 and cov >= 0.25
    moderate_protein = pid >= 0.65 and cov >= 0.15
    poor_protein = pid < 0.50 and cov < 0.15

    if strong_structure and strong_protein:
        candidate.fragment_class = "fragmented_high_confidence"
        candidate.classification_reason = "strong_structure_and_protein_support"
    elif moderate_structure and strong_protein:
        candidate.fragment_class = "fragmented_supported_partial"
        candidate.classification_reason = "partial_structure_but_strong_protein_support"
    elif moderate_structure and moderate_protein:
        candidate.fragment_class = "fragmented_moderate_confidence"
        candidate.classification_reason = "moderate_structure_and_protein_support"
    elif moderate_structure and poor_protein:
        candidate.fragment_class = "pseudogene_like_fragment"
        candidate.classification_reason = "structure_present_but_poor_protein_support"
    else:
        candidate.fragment_class = "fragmented_ambiguous"
        candidate.classification_reason = "mixed_or_insufficient_support"

    return candidate

def _max_consecutive_run(values: list[int]) -> int:
    if not values:
        return 0

    values = sorted(set(values))
    best = 1
    cur = 1

    for i in range(1, len(values)):
        if values[i] == values[i - 1] + 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1

    return best


def _safe_fraction(num: int, den: int) -> float:
    return num / den if den > 0 else 0.0

def classify_disrupted_projection_candidate(candidate):
    """
    Classify a fragmented/disrupted projection candidate.

    Expected fields on candidate:
      - n_target_seqids
      - exon_recovery_fraction
      - consecutive_exon_fraction
      - strand_consistent
      - exon_order_consistent
      - protein_support_class
    """
    prot = candidate.protein_support_class or "no_protein_evidence"
    rec = candidate.exon_recovery_fraction or 0.0
    cons = candidate.consecutive_exon_fraction or 0.0
    nseq = candidate.n_target_seqids or 0
    strand_ok = bool(candidate.strand_consistent)
    order_ok = bool(candidate.exon_order_consistent)

    # -----------------------------
    # Multi-scaffold split cases
    # -----------------------------
    if nseq >= 2:
        if strand_ok and order_ok and rec >= 0.90 and cons >= 0.90:
            if prot in {"strong", "moderate"}:
                candidate.fragment_class = "split_fragment_high_confidence"
                candidate.classification_reason = "multi_scaffold_coherent_with_good_protein_support"
            else:
                candidate.fragment_class = "split_fragment_structural"
                candidate.classification_reason = "multi_scaffold_coherent_but_protein_support_weak"
            return candidate

        if strand_ok and order_ok and rec >= 0.50:
            candidate.fragment_class = "split_fragment_ambiguous"
            candidate.classification_reason = "multi_scaffold_present_but_structure_incomplete"
            return candidate

        candidate.fragment_class = "ambiguous_partial"
        candidate.classification_reason = "multi_scaffold_mixed_support"
        return candidate

    # -----------------------------
    # Single-scaffold partial cases
    # -----------------------------
    if nseq == 1:
        if strand_ok and order_ok and rec >= 0.70 and cons >= 0.50:
            if prot in {"strong", "moderate"}:
                candidate.fragment_class = "partial_fragment_supported"
                candidate.classification_reason = "single_scaffold_coherent_partial_with_good_protein_support"
            elif prot == "weak":
                candidate.fragment_class = "ambiguous_partial"
                candidate.classification_reason = "single_scaffold_partial_structure_present_but_protein_support_weak"
            else:
                candidate.fragment_class = "ambiguous_partial"
                candidate.classification_reason = "single_scaffold_partial_mixed_support"
            return candidate

        candidate.fragment_class = "ambiguous_partial"
        candidate.classification_reason = "single_scaffold_partial_mixed_support"
        return candidate

    candidate.fragment_class = "ambiguous_partial"
    candidate.classification_reason = "insufficient_structure"
    return candidate
