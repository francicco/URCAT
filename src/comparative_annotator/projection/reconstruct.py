from __future__ import annotations

from collections import defaultdict

from comparative_annotator.models.projected_transcript import ProjectedTranscript


def _normalize_interval(start: int, end: int) -> tuple[int, int]:
    start_i = int(start)
    end_i = int(end)
    return (start_i, end_i) if start_i <= end_i else (end_i, start_i)


def _interval_length(start: int, end: int) -> int:
    start_n, end_n = _normalize_interval(start, end)
    return end_n - start_n + 1


def _detect_target_species(projected_exon_blocks, fallback: str | None = None) -> str:
    """
    Try to recover target species from HAL projection interval objects.

    This keeps backward compatibility with older callers that do not yet pass
    target_species explicitly.
    """
    for exon_hits in projected_exon_blocks:
        for iv in exon_hits:
            for attr in ("target_species", "species"):
                value = getattr(iv, attr, None)
                if value:
                    return value
    if fallback is not None:
        return fallback
    raise ValueError(
        "Could not determine target species in reconstruct_projected_transcripts(); "
        "pass target_species explicitly."
    )


def _group_hits_by_seqid_and_strand(projected_exon_blocks):
    """
    Flatten projected exon blocks and group them by (seqid, strand).

    Returns
    -------
    dict[(seqid, strand), list[dict]]
        Each row contains:
            source_exon_number
            target_start
            target_end
            chain_score
    """
    grouped = defaultdict(list)

    for exon_idx, exon_hits in enumerate(projected_exon_blocks, start=1):
        for iv in exon_hits:
            grouped[(iv.seqid, iv.strand)].append(
                {
                    "source_exon_number": exon_idx,
                    "target_start": int(iv.start),
                    "target_end": int(iv.end),
                    "chain_score": getattr(iv, "chain_score", None),
                }
            )

    return grouped


def _choose_best_hit_per_source_exon(rows: list[dict]) -> list[dict]:
    """
    For each source exon, keep the best target block.

    Current rule:
    - highest chain_score first
    - if tied or missing, prefer longer block
    - then earlier genomic coordinate for determinism
    """
    by_exon = defaultdict(list)
    for row in rows:
        by_exon[row["source_exon_number"]].append(row)

    chosen = []
    for exon_no in sorted(by_exon):
        best = sorted(
            by_exon[exon_no],
            key=lambda x: (
                -(x["chain_score"] if x["chain_score"] is not None else float("-inf")),
                -_interval_length(x["target_start"], x["target_end"]),
                _normalize_interval(x["target_start"], x["target_end"])[0],
                _normalize_interval(x["target_start"], x["target_end"])[1],
            ),
        )[0]
        chosen.append(best)

    return chosen


def _rows_to_exons(rows: list[dict], strand: str) -> list[tuple[int, int]]:
    """
    Convert chosen per-exon rows into ordered target exons.

    Genomic coordinates are always stored ascending as (start, end).
    Exon list is also sorted in genomic order, regardless of strand.
    """
    exons = [
        _normalize_interval(row["target_start"], row["target_end"])
        for row in rows
    ]
    exons.sort(key=lambda x: (x[0], x[1]))
    return exons


def _mean_chain_score(rows: list[dict]) -> float | None:
    vals = [float(r["chain_score"]) for r in rows if r.get("chain_score") is not None]
    if not vals:
        return None
    return sum(vals) / len(vals)


def reconstruct_projected_transcripts(
    seed_transcript,
    projected_exon_blocks,
    target_species: str | None = None,
    max_intron_gap: int | None = None,
) -> list[ProjectedTranscript]:
    """
    Reconstruct target-side projected transcripts from per-exon HAL projections.

    Parameters
    ----------
    seed_transcript
        Source transcript object. Expected fields:
            species, transcript_id, exons
    projected_exon_blocks
        List of lists of HAL projection intervals, one inner list per source exon.
    target_species
        Target species. Recommended to pass explicitly.
        If omitted, the function tries to infer it from the projection intervals.
    max_intron_gap
        Reserved for future use. Currently ignored.

    Returns
    -------
    list[ProjectedTranscript]
        One projected transcript per compatible (seqid, strand) group.
    """
    del max_intron_gap

    if not projected_exon_blocks:
        return []

    detected_target_species = _detect_target_species(
        projected_exon_blocks,
        fallback=target_species,
    )

    grouped = _group_hits_by_seqid_and_strand(projected_exon_blocks)
    projected_transcripts: list[ProjectedTranscript] = []

    n_source_exons = max(1, len(getattr(seed_transcript, "exons", []) or []))

    for (seqid, strand), rows in sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1])):
        chosen = _choose_best_hit_per_source_exon(rows)
        if not chosen:
            continue

        exons = _rows_to_exons(chosen, strand)
        if not exons:
            continue

        mean_chain = _mean_chain_score(chosen)
        mean_exon_recovery = len(chosen) / n_source_exons

        chain_orientation = "forward" if strand == getattr(seed_transcript, "strand", None) else "reverse"

        pt = ProjectedTranscript(
            species=detected_target_species,
            source_species=seed_transcript.species,
            source_transcript=seed_transcript.transcript_id,
            seqid=seqid,
            strand=strand,
            exons=exons,
            chain_orientation=chain_orientation,
            chain_score=mean_chain,
            mean_exon_recovery=mean_exon_recovery,
            support_score=mean_exon_recovery,
        )
        projected_transcripts.append(pt)

    return projected_transcripts
