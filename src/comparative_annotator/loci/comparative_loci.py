from __future__ import annotations

from comparative_annotator.models.comparative import ComparativeLocus
from comparative_annotator.models.projection import ProjectionInterval
from comparative_annotator.models.locus import SpeciesLocus

from comparative_annotator.projection.adjudication import choose_best_locus
from comparative_annotator.projection.matching import find_overlapping_loci

def overlaps(a_start, a_end, b_start, b_end):

    return a_start <= b_end and b_start <= a_end


def find_overlapping_loci(interval, species_loci):

    hits = []

    for locus in species_loci:

        if locus.seqid != interval.seqid:
            continue

        if locus.strand != interval.strand:
            continue

        if overlaps(interval.start, interval.end, locus.start, locus.end):
            hits.append(locus.locus_id)

    return hits


def build_comparative_locus(
    locus_id,
    seed_species,
    seed_transcript,
    projections,
    species_loci_dict,
):

    clocus = ComparativeLocus(
        locus_id=locus_id,
        seed_species=seed_species,
        seed_transcript=seed_transcript,
    )

    for proj in projections:

        loci = find_overlapping_loci(
            proj,
            species_loci_dict.get(proj.species, [])
        )

        for locus_id in loci:
            clocus.add_species_locus(proj.species, locus_id)

    return clocus
