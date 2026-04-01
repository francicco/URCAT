from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any

from comparative_annotator.models.projected_transcript import ProjectedTranscript


@dataclass
class _ProjectedBlock:
    source_exon_number: int
    target_seqid: str
    target_start: int
    target_end: int
    target_strand: str
    chain_score: float = 0.0


def _norm_interval(start: int, end: int) -> tuple[int, int]:
    return (start, end) if start <= end else (end, start)


def _interval_len(start: int, end: int) -> int:
    s, e = _norm_interval(start, end)
    return e - s + 1


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []

    ordered = sorted(_norm_interval(s, e) for s, e in intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _choose_best_block_per_exon(projected_exon_blocks: list[list[Any]]) -> list[_ProjectedBlock]:
    chosen: list[_ProjectedBlock] = []
    for exon_no, hits in enumerate(projected_exon_blocks, start=1):
        if not hits:
            continue
        best = sorted(
            hits,
            key=lambda x: (
                getattr(x, "chain_score", 0.0) or 0.0,
                _interval_len(int(x.start), int(x.end)),
            ),
            reverse=True,
        )[0]
        chosen.append(
            _ProjectedBlock(
                source_exon_number=exon_no,
                target_seqid=best.seqid,
                target_start=int(best.start),
                target_end=int(best.end),
                target_strand=best.strand,
                chain_score=float(getattr(best, "chain_score", 0.0) or 0.0),
            )
        )
    return chosen


def _group_blocks_into_transcripts(blocks: list[_ProjectedBlock]) -> list[list[_ProjectedBlock]]:
    by_seqid_strand: dict[tuple[str, str], list[_ProjectedBlock]] = defaultdict(list)
    for b in blocks:
        by_seqid_strand[(b.target_seqid, b.target_strand)].append(b)

    grouped: list[list[_ProjectedBlock]] = []
    for _, rows in by_seqid_strand.items():
        grouped.append(sorted(rows, key=lambda x: x.source_exon_number))
    return grouped


def reconstruct_projected_transcripts(
    seed_transcript,
    projected_exon_blocks: list[list[Any]],
    max_intron_gap: int | None = None,
) -> list[ProjectedTranscript]:
    """
    Reconstruct projected transcripts from exon-wise HAL projections.

    Current implementation chooses one best projection block per source exon,
    then groups compatible blocks by target seqid and strand.
    """
    if max_intron_gap is not None:
        raise NotImplementedError("max_intron_gap is not implemented in reconstruct_projected_transcripts")

    chosen_blocks = _choose_best_block_per_exon(projected_exon_blocks)
    if not chosen_blocks:
        return []

    txs: list[ProjectedTranscript] = []
    for rows in _group_blocks_into_transcripts(chosen_blocks):
        seqid = rows[0].target_seqid
        strand = rows[0].target_strand
        exon_intervals = _merge_intervals([(r.target_start, r.target_end) for r in rows])
        recovered_exons = [r.source_exon_number for r in rows]
        target_seqids = sorted({r.target_seqid for r in rows})
        total_chain_score = sum(r.chain_score for r in rows)
        mean_chain_score = total_chain_score / len(rows) if rows else 0.0

        txs.append(
            ProjectedTranscript(
                species=seed_transcript.species if False else getattr(seed_transcript, "target_species", None) or "",
                source_species=seed_transcript.species,
                source_transcript=seed_transcript.transcript_id,
                seqid=seqid,
                strand=strand,
                exons=exon_intervals,
                target_seqids=target_seqids,
                chain_orientation="forward" if strand == seed_transcript.strand else "reverse",
                total_chain_score=total_chain_score,
                mean_chain_score=mean_chain_score,
                exon_recovery_fraction=len(recovered_exons) / max(1, len(seed_transcript.exons)),
                n_source_exons=len(seed_transcript.exons),
                recovered_exon_numbers=recovered_exons,
                source_seqid=seed_transcript.seqid,
                source_strand=seed_transcript.strand,
            )
        )

    # Fill target species explicitly on the projected objects.
    for tx in txs:
        tx.species = getattr(seed_transcript, "_reconstruct_target_species", None) or tx.species

    return txs
