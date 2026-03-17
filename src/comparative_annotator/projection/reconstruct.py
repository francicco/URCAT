from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import List

from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.projection import ProjectionInterval


@dataclass
class IndexedProjection:
    source_exon_index: int
    interval: ProjectionInterval


def merge_projection_blocks(
    intervals: list[ProjectionInterval],
) -> list[ProjectionInterval]:
    """
    Merge fragmented projection blocks for a single source exon.

    Version 1:
    - merge all blocks that belong to the same species / seqid / strand / source transcript
    - do not impose a max gap threshold
    - preserve only the outer span
    """
    if not intervals:
        return []

    grouped: dict[tuple[str, str, str, str, str], list[ProjectionInterval]] = defaultdict(list)

    for iv in intervals:
        key = (
            iv.species,
            iv.seqid,
            iv.strand,
            iv.source_species,
            iv.source_transcript,
        )
        grouped[key].append(iv)

    merged: list[ProjectionInterval] = []

    for key, blocks in grouped.items():
        blocks = sorted(blocks, key=lambda x: (x.start, x.end))

        species, seqid, strand, source_species, source_transcript = key

        merged.append(
            ProjectionInterval(
                species=species,
                seqid=seqid,
                start=blocks[0].start,
                end=blocks[-1].end,
                strand=strand,
                source_species=source_species,
                source_transcript=source_transcript,
                identity=None,
                coverage=None,
            )
        )

    return merged


def block_midpoint(interval: ProjectionInterval) -> float:
    return (interval.start + interval.end) / 2.0


def chain_score(items: list[IndexedProjection], orientation: str) -> float:
    """
    Simple chain score:
    - reward number of unique source exons recovered
    - mildly penalize fragmentation
    """
    recovered = len(set(x.source_exon_index for x in items))
    fragmentation_penalty = 0.05 * max(0, len(items) - recovered)
    orientation_bonus = 0.01  # tiny constant so score is always numeric
    return recovered - fragmentation_penalty + orientation_bonus


def reconstruct_projected_transcripts(
    source_transcript,
    projected_exon_blocks: list[list[ProjectionInterval]],
) -> list[ProjectedTranscript]:
    """
    Reconstruct projected transcript candidates from exon-wise projections.

    New behavior:
    - merge fragmented blocks per source exon
    - group by target seqid and strand
    - allow both forward and reverse co-linear chains
    - keep only the best chain per target seqid/strand
    """
    indexed: list[IndexedProjection] = []

    for exon_idx, blocks in enumerate(projected_exon_blocks):
        merged_blocks = merge_projection_blocks(blocks)
        for block in merged_blocks:
            indexed.append(IndexedProjection(source_exon_index=exon_idx, interval=block))

    if not indexed:
        return []

    by_target: dict[tuple[str, str], list[IndexedProjection]] = defaultdict(list)
    for item in indexed:
        key = (item.interval.seqid, item.interval.strand)
        by_target[key].append(item)

    reconstructed: list[ProjectedTranscript] = []

    for (seqid, strand), items in by_target.items():
        forward_items = sorted(items, key=lambda x: (block_midpoint(x.interval), x.source_exon_index))
        reverse_items = sorted(items, key=lambda x: (-block_midpoint(x.interval), x.source_exon_index))

        forward_chain = _build_colinear_chain(
            forward_items,
            source_species=source_transcript.species,
            source_transcript_id=source_transcript.transcript_id,
            target_species=items[0].interval.species,
            seqid=seqid,
            strand=strand,
            orientation="forward",
        )

        reverse_chain = _build_colinear_chain(
            reverse_items,
            source_species=source_transcript.species,
            source_transcript_id=source_transcript.transcript_id,
            target_species=items[0].interval.species,
            seqid=seqid,
            strand=strand,
            orientation="reverse",
        )

        candidates = [c for c in [forward_chain, reverse_chain] if c is not None]
        if not candidates:
            continue

        # Keep only the best candidate for this target seqid/strand.
        # Prefer more recovered source exons; break ties with chain_score.
        best = max(
            candidates,
            key=lambda c: (
                c.coverage if c.coverage is not None else 0,
                c.chain_score if c.chain_score is not None else 0.0,
            ),
        )

        reconstructed.append(best)

    reconstructed.sort(key=lambda x: (x.chain_score or 0), reverse=True)
    return reconstructed


def _build_colinear_chain(
    items: list[IndexedProjection],
    source_species: str,
    source_transcript_id: str,
    target_species: str,
    seqid: str,
    strand: str,
    orientation: str,
) -> ProjectedTranscript | None:
    """
    Build one co-linear chain assuming target coordinates are already ordered
    in either forward or reverse genomic direction.

    Requirement:
    - source exon indices must increase monotonically
    """
    if not items:
        return None

    selected: list[IndexedProjection] = []
    seen_exons: set[int] = set()
    last_exon_idx = -1

    for item in items:
        exon_idx = item.source_exon_index

        if exon_idx < last_exon_idx:
            continue

        # keep only first compatible projection per source exon in this simple version
        if exon_idx in seen_exons:
            continue

        selected.append(item)
        seen_exons.add(exon_idx)
        last_exon_idx = exon_idx

    if not selected:
        return None

    pt = ProjectedTranscript(
        species=target_species,
        seqid=seqid,
        strand=strand,
        source_species=source_species,
        source_transcript=source_transcript_id,
    )

    for item in selected:
        pt.add_exon(
            item.interval.start,
            item.interval.end,
            source_exon_index=item.source_exon_index,
        )

    # Sort exons in genomic order for storage/printing
    pt.exons.sort()
    pt.coverage = len(set(pt.source_exon_indices))
    pt.chain_orientation = orientation
    pt.fragmentation_count = len(items) - len(set(x.source_exon_index for x in items))
    pt.chain_score = chain_score(selected, orientation=orientation)

    return pt
