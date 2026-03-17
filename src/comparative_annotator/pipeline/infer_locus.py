from __future__ import annotations

from comparative_annotator.loci.comparative_builder import (
    build_comparative_locus_from_projection,
)
from comparative_annotator.projection.reconstruct import (
    reconstruct_projected_transcripts,
)


def infer_comparative_locus(
    seed_transcript,
    target_species,
    hal_adapter,
    species_loci,
    transcripts_by_species=None,
):
    projected_exon_blocks = []

    for exon_start, exon_end in seed_transcript.exons:
        intervals = hal_adapter.project_interval(
            source_species=seed_transcript.species,
            target_species=target_species,
            seqid=seed_transcript.seqid,
            start=exon_start,
            end=exon_end,
            strand=seed_transcript.strand,
            source_transcript=seed_transcript.transcript_id,
        )
        projected_exon_blocks.append(intervals)

    projected_transcripts = reconstruct_projected_transcripts(
        seed_transcript,
        projected_exon_blocks,
    )

    clocus = build_comparative_locus_from_projection(
        seed_species=seed_transcript.species,
        seed_transcript=seed_transcript.transcript_id,
        projected_transcripts=projected_transcripts,
        species_loci=species_loci,
        source_transcript=seed_transcript,
        transcripts_by_species=transcripts_by_species,
    )

    return clocus
