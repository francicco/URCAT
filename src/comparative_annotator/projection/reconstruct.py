from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, List

from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.projection import ProjectionInterval


@dataclass
class IndexedProjection:
    source_exon_index: int
    interval: ProjectionInterval


def reconstruct_projected_transcripts(
    source_transcript,
    projected_exon_blocks: list[list[ProjectionInterval]],
    max_exon_gap: int = 1_000_000,
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
    max_exon_gap
        Very permissive distance threshold used to split implausibly distant chains.

    Returns
    -------
    list[ProjectedTranscript]
        One or more projected transcript candidates.
    """
    indexed: list[IndexedProjection] = []

    for exon_idx, blocks in enumerate(projected_exon_blocks):
        for block in blocks:
            indexed.append(IndexedProjection(source_exon_index=exon_idx, interval=block))

    if not indexed:
        return []

    # Group by target seqid and strand first
    by_target: dict[tuple[str, str], list[IndexedProjection]] = defaultdict(list)
    for item in indexed:
        key = (item.interval.seqid, item.interval.strand)
        by_target[key].append(item)

    reconstructed: list[ProjectedTranscript] = []

    for (seqid, strand), items in by_target.items():
        # Sort by genomic position, then by source exon order
        items.sort(key=lambda x: (x.interval.start, x.interval.end, x.source_exon_index))

        chains = _split_into_compatible_chains(
            items,
            source_species=source_transcript.species,
            source_transcript_id=source_transcript.transcript_id,
            target_species=items[0].interval.species,
            seqid=seqid,
            strand=strand,
            max_exon_gap=max_exon_gap,
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
    max_exon_gap: int,
) -> list[ProjectedTranscript]:
    """
    Split projected exon blocks into one or more compatible chains.

    Version 1 logic:
    - projected exon indices should not go backwards within a chain
    - genomic coordinates should move forward
    - huge genomic jumps create a new chain
    """
    if not items:
        return []

    chains: list[list[IndexedProjection]] = []
    current_chain: list[IndexedProjection] = [items[0]]

    last_exon_idx = items[0].source_exon_index
    last_start = items[0].interval.start
    last_end = items[0].interval.end

    for item in items[1:]:
        exon_idx = item.source_exon_index
        start = item.interval.start
        end = item.interval.end

        same_or_forward_source_order = exon_idx >= last_exon_idx
        forward_genomic_order = start >= last_start
        gap_ok = (start - last_end) <= max_exon_gap

        if same_or_forward_source_order and forward_genomic_order and gap_ok:
            current_chain.append(item)
        else:
            chains.append(current_chain)
            current_chain = [item]

        last_exon_idx = exon_idx
        last_start = start
        last_end = end

    if current_chain:
        chains.append(current_chain)

    output: list[ProjectedTranscript] = []

    for chain_idx, chain in enumerate(chains, start=1):
        pt = ProjectedTranscript(
            species=target_species,
            seqid=seqid,
            strand=strand,
            source_species=source_species,
            source_transcript=source_transcript_id,
        )

        # Keep exon order as recovered in this chain
        for item in chain:
            pt.add_exon(item.interval.start, item.interval.end)

        # Store useful metadata
        pt.exons.sort()
        recovered_source_exons = [item.source_exon_index for item in chain]
        pt.coverage = len(set(recovered_source_exons))
        pt.identity = None

        output.append(pt)

    return output
