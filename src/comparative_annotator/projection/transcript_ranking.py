from __future__ import annotations

from comparative_annotator.models.transcript_match import ProjectionTranscriptScore


def span_overlap_fraction(projected, target_transcript) -> float:
    if not projected.exons or not target_transcript.exons:
        return 0.0

    p_start = min(s for s, e in projected.exons)
    p_end = max(e for s, e in projected.exons)
    t_start = min(s for s, e in target_transcript.exons)
    t_end = max(e for s, e in target_transcript.exons)

    overlap = max(0, min(p_end, t_end) - max(p_start, t_start) + 1)
    p_len = p_end - p_start + 1

    if p_len <= 0:
        return 0.0

    return overlap / p_len


def exon_count_similarity(projected, target_transcript) -> float:
    p = max(1, projected.exon_count)
    t = max(1, target_transcript.exon_count)
    return min(p, t) / max(p, t)


def orientation_consistency(projected, target_transcript) -> float:
    return 1.0 if projected.strand == target_transcript.strand else 0.7


def intron_support(projected, target_transcript) -> float:
    """
    Very light version 1.
    Reward matching intron counts, without requiring exact coordinates yet.
    """
    p_introns = max(0, projected.exon_count - 1)
    t_introns = getattr(target_transcript, "intron_count", max(0, target_transcript.exon_count - 1))

    if p_introns == 0 and t_introns == 0:
        return 1.0

    denom = max(1, max(p_introns, t_introns))
    return min(p_introns, t_introns) / denom


def exon_recovery(projected, source_transcript) -> float:
    if not source_transcript.exons:
        return 0.0
    recovered = len(projected.source_exon_indices or [])
    return recovered / len(source_transcript.exons)


def score_projected_transcript_against_transcript(
    projected,
    source_transcript,
    target_transcript,
    locus_id: str | None = None,
) -> ProjectionTranscriptScore:
    er = exon_recovery(projected, source_transcript)
    ecs = exon_count_similarity(projected, target_transcript)
    so = span_overlap_fraction(projected, target_transcript)
    oc = orientation_consistency(projected, target_transcript)
    ins = intron_support(projected, target_transcript)

    total = (
        0.30 * er +
        0.20 * ecs +
        0.25 * so +
        0.10 * oc +
        0.15 * ins
    )

    return ProjectionTranscriptScore(
        species=projected.species,
        projected_id=projected.source_transcript,
        transcript_id=target_transcript.transcript_id,
        locus_id=locus_id,
        exon_recovery=er,
        exon_count_similarity=ecs,
        span_overlap=so,
        orientation_consistency=oc,
        intron_support=ins,
        total_score=total,
    )


def rank_transcripts_within_locus(
    projected,
    source_transcript,
    locus,
    target_transcripts_by_id,
):
    scores = []

    for tx_id in locus.transcripts:
        if tx_id not in target_transcripts_by_id:
            continue

        tx = target_transcripts_by_id[tx_id]
        score = score_projected_transcript_against_transcript(
            projected=projected,
            source_transcript=source_transcript,
            target_transcript=tx,
            locus_id=locus.locus_id,
        )
        scores.append(score)

    scores.sort(key=lambda x: x.total_score, reverse=True)
    return scores


def choose_best_transcript_within_locus(
    projected,
    source_transcript,
    locus,
    target_transcripts_by_id,
):
    ranked = rank_transcripts_within_locus(
        projected=projected,
        source_transcript=source_transcript,
        locus=locus,
        target_transcripts_by_id=target_transcripts_by_id,
    )

    if not ranked:
        return None, []

    return ranked[0], ranked[1:]
