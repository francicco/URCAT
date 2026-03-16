from __future__ import annotations

from pprint import pprint

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
):
    """
    Infer a comparative locus for one seed transcript against one target species.

    Steps:
    1. project each exon with HAL
    2. reconstruct projected transcript candidates
    3. match/adjudicate against target species loci
    4. return ComparativeLocus
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
        print(f"Projected exon {exon_start}-{exon_end} to {target_species}:")
        pprint(intervals)
        projected_exon_blocks.append(intervals)

    projected_transcripts = reconstruct_projected_transcripts(
        seed_transcript,
        projected_exon_blocks,
    )

    print(f"Reconstructed projected transcripts for {target_species}:")
    pprint(projected_transcripts)

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
    )

    return clocus
