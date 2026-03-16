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

    This is deliberately permissive. Later versions can add
    stricter fragmentation logic or CESAR-like refinement.
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


def reconstruct_projected_transcripts(
    source_transcript,
    projected_exon_blocks: list[list[ProjectionInterval]],
) -> list[ProjectedTranscript]:
    """
    Reconstruct projected transcript candidates from exon-wise projections.

    Parameters
    ----------
    source_transcript
        A CandidateTranscript-like object with transcript_id, species, strand, exons.
    projected_exon_blocks
        A list with one entry per source exon.
        Each entry is a list of ProjectionInterval objects produced by projecting that exon.

    Returns
    -------
    list[ProjectedTranscript]
        One or more projected transcript candidates.
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
        items.sort(key=lambda x: (x.interval.start, x.interval.end, x.source_exon_index))

        chains = _split_into_compatible_chains(
            items,
            source_species=source_transcript.species,
            source_transcript_id=source_transcript.transcript_id,
            target_species=items[0].interval.species,
            seqid=seqid,
            strand=strand,
        )

        reconstructed.extend(chains)

    return reconstructed


def _split_into_compatible_chains(
    items: list[IndexedProjection],
    source_species: str,
    source_transcript_id: str,
    target_species: str,
    seqid: str,
    strand: str,
) -> list[ProjectedTranscript]:
    """
    Split projected exon blocks into one or more compatible chains.

    Version 1 logic:
    - source exon indices should not go backwards within a chain
    - genomic coordinates should move forward
    - no explicit max genomic gap threshold is imposed
    """
    if not items:
        return []

    chains: list[list[IndexedProjection]] = []
    current_chain: list[IndexedProjection] = [items[0]]

    last_exon_idx = items[0].source_exon_index
    last_start = items[0].interval.start

    for item in items[1:]:
        exon_idx = item.source_exon_index
        start = item.interval.start

        same_or_forward_source_order = exon_idx >= last_exon_idx
        forward_genomic_order = start >= last_start

        if same_or_forward_source_order and forward_genomic_order:
            current_chain.append(item)
        else:
            chains.append(current_chain)
            current_chain = [item]

        last_exon_idx = exon_idx
        last_start = start

    if current_chain:
        chains.append(current_chain)

    output: list[ProjectedTranscript] = []

    for chain in chains:
        pt = ProjectedTranscript(
            species=target_species,
            seqid=seqid,
            strand=strand,
            source_species=source_species,
            source_transcript=source_transcript_id,
        )

        recovered_source_exons = []

        for item in chain:
            pt.add_exon(item.interval.start, item.interval.end)
            recovered_source_exons.append(item.source_exon_index)

        pt.exons.sort()
        pt.coverage = len(set(recovered_source_exons))

        output.append(pt)

    return output
