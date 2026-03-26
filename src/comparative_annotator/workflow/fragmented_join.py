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
    candidate: FragmentedComparativeLocus,
    *,
    diamond_hits_for_source: dict | None,
    source_protein_length: int | None = None,
) -> FragmentedComparativeLocus:
    """
    diamond_hits_for_source:
      dict of {(qseqid, sseqid): {"pid": ..., "aln_len": ..., "bitscore": ..., "evalue": ...}}
      typically one source-vs-target DIAMOND result table.

    We summarize support over all hits for the source transcript.
    """
    if not diamond_hits_for_source:
        candidate.protein_support_class = "no_protein_evidence"
        return candidate

    source_tx = candidate.source_transcript
    matched = []
    for (qseqid, sseqid), hit in diamond_hits_for_source.items():
        if qseqid == source_tx:
            matched.append(hit)

    if not matched:
        candidate.protein_support_class = "no_protein_evidence"
        return candidate

    matched = sorted(matched, key=lambda x: x.get("bitscore", 0.0), reverse=True)

    candidate.aa_identity_mean = _mean([float(h["pid"]) for h in matched]) / 100.0
    candidate.aa_bitscore = float(matched[0]["bitscore"])

    if source_protein_length and source_protein_length > 0:
        best_aln_len = max(int(h["aln_len"]) for h in matched)
        candidate.aa_coverage_fraction = best_aln_len / source_protein_length
    else:
        candidate.aa_coverage_fraction = None

    pid = candidate.aa_identity_mean or 0.0
    cov = candidate.aa_coverage_fraction if candidate.aa_coverage_fraction is not None else 0.0

    if pid >= 0.80 and cov >= 0.30:
        candidate.protein_support_class = "strong"
    elif pid >= 0.65 and cov >= 0.20:
        candidate.protein_support_class = "moderate"
    elif pid >= 0.50:
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

def classify_disrupted_projection_candidate(
    candidate: FragmentedComparativeLocus,
) -> FragmentedComparativeLocus:
    """
    Classify a disrupted projection into one of:
      - split_fragment_high_confidence
      - split_fragment_moderate_confidence
      - partial_ortholog_high_confidence
      - partial_ortholog_moderate_confidence
      - spurious_partial
      - ambiguous_partial

    Logic:
      1. Separate multi-scaffold split cases from single-scaffold partial cases.
      2. For single-scaffold partials, distinguish coherent partial orthologs
         from noisy/spurious partial mappings.
    """

    recovered = sorted(set(candidate.recovered_exon_numbers))
    n_recovered = len(recovered)
    n_source = max(1, candidate.n_source_exons)

    max_consecutive = _max_consecutive_run(recovered)
    consecutive_fraction = _safe_fraction(max_consecutive, n_source)

    # store if model has these attributes
    if hasattr(candidate, "max_consecutive_recovered_exons"):
        candidate.max_consecutive_recovered_exons = max_consecutive
    if hasattr(candidate, "consecutive_exon_fraction"):
        candidate.consecutive_exon_fraction = consecutive_fraction

    pid = candidate.aa_identity_mean or 0.0
    cov = candidate.aa_coverage_fraction if candidate.aa_coverage_fraction is not None else 0.0
    n_seqids = len(candidate.target_seqids)

    strong_protein = pid >= 0.80 and cov >= 0.25
    moderate_protein = pid >= 0.65 and cov >= 0.15
    weak_protein = pid >= 0.50
    poor_protein = pid < 0.50 and cov < 0.15

    strong_structure = (
        candidate.exon_recovery_fraction >= 0.40
        and n_recovered >= 2
        and candidate.strand_consistent
        and candidate.exon_order_consistent
    )

    coherent_partial_structure = (
        candidate.exon_recovery_fraction >= 0.20
        and n_recovered >= 2
        and candidate.strand_consistent
        and candidate.exon_order_consistent
        and max_consecutive >= 2
    )

    very_weak_structure = (
        n_recovered <= 1
        or not candidate.strand_consistent
        or not candidate.exon_order_consistent
    )

    # ------------------------------------------------------------------
    # A. Multi-scaffold split cases
    # ------------------------------------------------------------------
    if n_seqids >= 2:
        if strong_structure and strong_protein:
            candidate.fragment_class = "split_fragment_high_confidence"
            candidate.classification_reason = "multi_scaffold_coherent_with_strong_protein_support"
            return candidate

        if coherent_partial_structure and moderate_protein:
            candidate.fragment_class = "split_fragment_moderate_confidence"
            candidate.classification_reason = "multi_scaffold_coherent_with_moderate_protein_support"
            return candidate

        if poor_protein and very_weak_structure:
            candidate.fragment_class = "spurious_partial"
            candidate.classification_reason = "multi_scaffold_but_weak_structure_and_protein_support"
            return candidate

        candidate.fragment_class = "ambiguous_partial"
        candidate.classification_reason = "multi_scaffold_mixed_support"
        return candidate

    # ------------------------------------------------------------------
    # B. Single-scaffold partial cases
    # ------------------------------------------------------------------
    # Here the question is:
    # is this a real truncated/incomplete ortholog or just a spurious patch?
    if n_seqids == 1:
        if coherent_partial_structure and strong_protein and consecutive_fraction >= 0.15:
            candidate.fragment_class = "partial_ortholog_high_confidence"
            candidate.classification_reason = "single_scaffold_coherent_partial_with_strong_protein_support"
            return candidate

        if coherent_partial_structure and moderate_protein and max_consecutive >= 2:
            candidate.fragment_class = "partial_ortholog_moderate_confidence"
            candidate.classification_reason = "single_scaffold_coherent_partial_with_moderate_protein_support"
            return candidate

        if poor_protein and (very_weak_structure or consecutive_fraction < 0.10):
            candidate.fragment_class = "spurious_partial"
            candidate.classification_reason = "single_scaffold_partial_with_weak_structure_and_protein_support"
            return candidate

        if weak_protein and coherent_partial_structure:
            candidate.fragment_class = "ambiguous_partial"
            candidate.classification_reason = "single_scaffold_partial_structure_present_but_protein_support_weak"
            return candidate

        candidate.fragment_class = "ambiguous_partial"
        candidate.classification_reason = "single_scaffold_partial_mixed_support"
        return candidate

    # ------------------------------------------------------------------
    # C. Fallback
    # ------------------------------------------------------------------
    candidate.fragment_class = "ambiguous_partial"
    candidate.classification_reason = "unclassified_disrupted_projection"
    return candidate
