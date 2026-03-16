from __future__ import annotations

from comparative_annotator.models.reciprocal import ReciprocalProjectionResult
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus


def validate_reciprocal_projection(
    seed_transcript,
    source_species_locus_id: str,
    target_species: str,
    hal_adapter,
    species_loci,
    representative_transcript_getter,
):
    """
    Validate reciprocal projection consistency.

    Parameters
    ----------
    seed_transcript
        Transcript in the source species.
    source_species_locus_id
        Locus containing the seed transcript in the source species.
    target_species
        Species to project into.
    hal_adapter
        HALAdapter instance.
    species_loci
        Dict: species -> list[SpeciesLocus]
    representative_transcript_getter
        Callable taking (species, locus_id) and returning a representative transcript object.

    Returns
    -------
    ReciprocalProjectionResult
    """
    forward = infer_comparative_locus(
        seed_transcript=seed_transcript,
        target_species=target_species,
        hal_adapter=hal_adapter,
        species_loci=species_loci,
    )

    forward_primary = forward.primary.get(target_species)

    if forward_primary is None:
        return ReciprocalProjectionResult(
            source_species=seed_transcript.species,
            target_species=target_species,
            seed_transcript=seed_transcript.transcript_id,
            forward_primary_locus=None,
            reverse_primary_locus=None,
            reciprocal_supported=False,
            classification="missing_forward_match",
        )

    target_representative = representative_transcript_getter(target_species, forward_primary)

    reverse = infer_comparative_locus(
        seed_transcript=target_representative,
        target_species=seed_transcript.species,
        hal_adapter=hal_adapter,
        species_loci=species_loci,
    )

    reverse_primary = reverse.primary.get(seed_transcript.species)

    if reverse_primary is None:
        return ReciprocalProjectionResult(
            source_species=seed_transcript.species,
            target_species=target_species,
            seed_transcript=seed_transcript.transcript_id,
            forward_primary_locus=forward_primary,
            reverse_primary_locus=None,
            reciprocal_supported=False,
            classification="missing_reverse_match",
        )

    reciprocal_supported = reverse_primary == source_species_locus_id

    return ReciprocalProjectionResult(
        source_species=seed_transcript.species,
        target_species=target_species,
        seed_transcript=seed_transcript.transcript_id,
        forward_primary_locus=forward_primary,
        reverse_primary_locus=reverse_primary,
        reciprocal_supported=reciprocal_supported,
        classification=(
            "reciprocal_ortholog"
            if reciprocal_supported
            else "nonreciprocal_candidate"
        ),
    )
