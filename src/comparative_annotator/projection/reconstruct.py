from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ProjectedTranscript:
    source_species: str
    target_species: str
    source_transcript: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    source_exon_numbers: list[int] = field(default_factory=list)
    source_blocks: list[dict[str, Any]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(s for s, _ in self.exons) if self.exons else 0

    @property
    def end(self) -> int:
        return max(e for _, e in self.exons) if self.exons else 0


@dataclass
class _ProjectedBlock:
    source_exon_number: int
    target_seqid: str
    target_start: int
    target_end: int
    target_strand: str
    chain_score: float | None = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _norm_interval(start: int, end: int) -> tuple[int, int]:
    return (start, end) if start <= end else (end, start)


def _interval_len(start: int, end: int) -> int:
    start, end = _norm_interval(start, end)
    return end - start + 1


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


def _coerce_projection_block(iv: Any, source_exon_number: int) -> _ProjectedBlock:
    return _ProjectedBlock(
        source_exon_number=int(source_exon_number),
        target_seqid=iv.seqid,
        target_start=int(iv.start),
        target_end=int(iv.end),
        target_strand=iv.strand,
        chain_score=getattr(iv, "chain_score", None),
    )


def _flatten_projected_exon_blocks(projected_exon_blocks: list[list[Any]]) -> list[_ProjectedBlock]:
    out: list[_ProjectedBlock] = []
    for exon_idx, hits in enumerate(projected_exon_blocks, start=1):
        for iv in hits:
            out.append(_coerce_projection_block(iv, exon_idx))
    return out


def _dominant_strand(rows: list[_ProjectedBlock], fallback: str) -> str:
    support: dict[str, int] = {}
    for r in rows:
        support[r.target_strand] = support.get(r.target_strand, 0) + 1
    if not support:
        return fallback
    return sorted(support.items(), key=lambda x: (-x[1], x[0]))[0][0]


def _group_blocks_by_seqid_and_strand(rows: list[_ProjectedBlock]) -> dict[tuple[str, str], list[_ProjectedBlock]]:
    grouped: dict[tuple[str, str], list[_ProjectedBlock]] = {}
    for r in rows:
        grouped.setdefault((r.target_seqid, r.target_strand), []).append(r)
    return grouped


def _pick_best_row_per_source_exon(rows: list[_ProjectedBlock]) -> list[_ProjectedBlock]:
    by_exon: dict[int, list[_ProjectedBlock]] = {}
    for r in rows:
        by_exon.setdefault(r.source_exon_number, []).append(r)

    chosen: list[_ProjectedBlock] = []
    for exon_no in sorted(by_exon):
        exon_rows = sorted(
            by_exon[exon_no],
            key=lambda x: ((x.chain_score or 0.0), _interval_len(x.target_start, x.target_end)),
            reverse=True,
        )
        chosen.append(exon_rows[0])
    return chosen


def _coalesce_same_exon_fragments(rows: list[_ProjectedBlock]) -> list[_ProjectedBlock]:
    by_key: dict[tuple[int, str, str], list[tuple[int, int]]] = {}
    chain_scores: dict[tuple[int, str, str], float | None] = {}

    for r in rows:
        key = (r.source_exon_number, r.target_seqid, r.target_strand)
        by_key.setdefault(key, []).append((r.target_start, r.target_end))
        best_score = chain_scores.get(key)
        cur_score = r.chain_score
        if best_score is None or ((cur_score or 0.0) > (best_score or 0.0)):
            chain_scores[key] = cur_score

    out: list[_ProjectedBlock] = []
    for (exon_no, seqid, strand), intervals in by_key.items():
        for start, end in _merge_intervals(intervals):
            out.append(
                _ProjectedBlock(
                    source_exon_number=exon_no,
                    target_seqid=seqid,
                    target_start=start,
                    target_end=end,
                    target_strand=strand,
                    chain_score=chain_scores[(exon_no, seqid, strand)],
                )
            )
    return sorted(out, key=lambda x: (x.source_exon_number, x.target_start, x.target_end))


def _rows_to_projected_transcript(
    rows: list[_ProjectedBlock],
    seed_transcript: Any,
    target_species: str,
) -> ProjectedTranscript | None:
    if not rows:
        return None

    rows = _coalesce_same_exon_fragments(rows)
    seqid = rows[0].target_seqid
    strand = rows[0].target_strand or getattr(seed_transcript, "strand", "+")

    exons = [_norm_interval(r.target_start, r.target_end) for r in rows]
    source_exon_numbers = [r.source_exon_number for r in rows]
    source_blocks = [
        {
            "source_exon_number": r.source_exon_number,
            "target_seqid": r.target_seqid,
            "target_start": r.target_start,
            "target_end": r.target_end,
            "target_strand": r.target_strand,
            "chain_score": r.chain_score,
        }
        for r in rows
    ]

    return ProjectedTranscript(
        source_species=getattr(seed_transcript, "species", "unknown"),
        target_species=target_species,
        source_transcript=getattr(seed_transcript, "transcript_id", getattr(seed_transcript, "source_transcript", "unknown")),
        seqid=seqid,
        strand=strand,
        exons=exons,
        source_exon_numbers=source_exon_numbers,
        source_blocks=source_blocks,
    )


def _target_order_consistent(rows: list[_ProjectedBlock]) -> bool:
    if not rows:
        return True
    strand = rows[0].target_strand
    rows = sorted(rows, key=lambda x: x.source_exon_number)
    starts = [_norm_interval(r.target_start, r.target_end)[0] for r in rows]
    if strand == "-":
        return starts == sorted(starts, reverse=True)
    return starts == sorted(starts)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def reconstruct_projected_transcripts(
    seed_transcript: Any,
    projected_exon_blocks: list[list[Any]],
    max_intron_gap: int | None = None,
) -> list[ProjectedTranscript]:
    """
    Reconstruct conservative projected transcript models from flattened HAL exon projections.

    Current behavior is intentionally strict:
    - blocks are separated by target seqid and strand
    - within each (seqid, strand) group, only the best block per source exon is kept
    - a group is emitted only if exon order on the target is consistent
    - each emitted transcript is a projected transcript candidate, not a fully validated orthology call

    Parameters
    ----------
    seed_transcript:
        Source transcript object. Must expose at least `species`, `transcript_id`, and `strand`.
    projected_exon_blocks:
        One list of HAL intervals per source exon, in source exon order.
    max_intron_gap:
        Reserved for later use. Currently ignored.
    """
    del max_intron_gap

    flat_rows = _flatten_projected_exon_blocks(projected_exon_blocks)
    if not flat_rows:
        return []

    grouped = _group_blocks_by_seqid_and_strand(flat_rows)
    transcripts: list[ProjectedTranscript] = []

    for (target_seqid, target_strand), rows in sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1])):
        chosen = _pick_best_row_per_source_exon(rows)
        if not chosen:
            continue

        if not _target_order_consistent(chosen):
            # keep the candidate only if at least two exons still support a usable chain after monotonic filtering
            monotonic: list[_ProjectedBlock] = []
            chosen_sorted = sorted(chosen, key=lambda x: x.source_exon_number)
            last_pos: int | None = None
            for r in chosen_sorted:
                pos = _norm_interval(r.target_start, r.target_end)[0]
                if last_pos is None:
                    monotonic.append(r)
                    last_pos = pos
                    continue
                if target_strand == "-":
                    if pos <= last_pos:
                        monotonic.append(r)
                        last_pos = pos
                else:
                    if pos >= last_pos:
                        monotonic.append(r)
                        last_pos = pos
            chosen = monotonic

        if not chosen:
            continue

        tx = _rows_to_projected_transcript(
            rows=chosen,
            seed_transcript=seed_transcript,
            target_species=getattr(rows[0], "target_species", None) or "unknown",
        )
        if tx is not None:
            tx.seqid = target_seqid
            tx.strand = _dominant_strand(chosen, getattr(seed_transcript, "strand", "+"))
            transcripts.append(tx)

    transcripts.sort(key=lambda tx: (tx.seqid, tx.start, tx.end, tx.source_transcript))
    return transcripts
