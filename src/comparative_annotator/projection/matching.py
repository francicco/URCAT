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
    Return all loci overlapping a projected transcript candidate.
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


def locus_overlap_fraction(
    projected: ProjectedTranscript,
    locus: SpeciesLocus,
):
    """
    Fraction of the projected transcript span overlapping a locus.
    """

    overlap_start = max(projected.start, locus.start)
    overlap_end = min(projected.end, locus.end)

    overlap = max(0, overlap_end - overlap_start)

    proj_len = projected.end - projected.start

    if proj_len == 0:
        return 0

    return overlap / proj_len
