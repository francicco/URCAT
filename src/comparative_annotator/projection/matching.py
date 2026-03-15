from __future__ import annotations

from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus


def overlaps(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


def find_overlapping_species_loci(
    projected: ProjectedTranscript,
    species_loci: list[SpeciesLocus],
):
    """
    Find loci overlapping a projected transcript candidate.
    """

    hits = []

    for locus in species_loci:

        if locus.seqid != projected.seqid:
            continue

        if locus.strand != projected.strand:
            continue

        if overlaps(projected.start, projected.end, locus.start, locus.end):
            hits.append(locus)

    return hits
