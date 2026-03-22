from __future__ import annotations

from collections import defaultdict


def _get(block, key, default=None):
    if isinstance(block, dict):
        return block.get(key, default)
    return getattr(block, key, default)


def summarize_projected_blocks(
    source_species,
    target_species,
    source_transcripts,
    blocks,
):
    """
    Summarise projected exon blocks for one source locus/transcript set into a
    fragmentation-aware record.

    Accepts blocks either as:
      - dicts with keys like 'target_seqid', 'target_start', ...
      - objects with attributes like .target_seqid, .target_start, ...

    Returns a dict suitable for JSON serialisation / table writing.
    """

    by_seqid = defaultdict(list)
    for b in blocks:
        target_seqid = _get(b, "target_seqid")
        if target_seqid is None:
            continue
        by_seqid[target_seqid].append(b)

    target_seqids = sorted(by_seqid.keys())
    n_target_seqids = len(target_seqids)

    all_intervals = []
    for seqid in target_seqids:
        seq_blocks = sorted(
            by_seqid[seqid],
            key=lambda x: (_get(x, "target_start", 0), _get(x, "target_end", 0)),
        )
        for b in seq_blocks:
            all_intervals.append(
                {
                    "target_seqid": _get(b, "target_seqid"),
                    "target_start": _get(b, "target_start"),
                    "target_end": _get(b, "target_end"),
                    "target_strand": _get(b, "target_strand"),
                    "source_exon_number": _get(b, "source_exon_number"),
                    "chain_score": _get(b, "chain_score"),
                }
            )

    is_fragmented = n_target_seqids > 1

    return {
        "source_species": source_species,
        "target_species": target_species,
        "source_transcripts": sorted(source_transcripts),
        "n_source_transcripts": len(source_transcripts),
        "n_projected_blocks": len(blocks),
        "n_target_seqids": n_target_seqids,
        "target_seqids": target_seqids,
        "is_fragmented": is_fragmented,
        "projected_intervals": all_intervals,
    }
