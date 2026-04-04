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
    """
    Project one seed transcript into one target species and build a ComparativeLocus.

    Parameters
    ----------
    seed_transcript
        Source transcript object. Expected to provide at least:
        - species
        - transcript_id
        - seqid
        - strand
        - exons
    target_species
        Species into which the seed transcript is projected.
    hal_adapter
        HAL adapter exposing `project_interval(...)`.
    species_loci
        Mapping of species -> locus objects, used by the comparative builder.
    transcripts_by_species
        Optional mapping of species -> transcript objects, forwarded to the
        comparative builder.

    Returns
    -------
    ComparativeLocus
        Locus-level interpretation of the projection in the target species.
    """
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
        seed_transcript=seed_transcript,
        projected_exon_blocks=projected_exon_blocks,
        target_species=target_species,
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
