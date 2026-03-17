from __future__ import annotations

from comparative_annotator.models.comparative import ComparativeLocus
from comparative_annotator.projection.adjudication import choose_best_locus
from comparative_annotator.projection.matching import (
    find_overlapping_species_loci,
    find_overlapping_species_loci_any_strand,
)


def build_comparative_locus_from_projection(
    seed_species,
    seed_transcript,
    projected_transcripts,
    species_loci,
):
    clocus = ComparativeLocus(
        locus_id=f"{seed_species}:{seed_transcript}",
        seed_species=seed_species,
        seed_transcript=seed_transcript,
    )

    for proj in projected_transcripts:
        loci_for_species = species_loci.get(proj.species, [])

        # 1. same-strand overlaps: normal adjudication
        same_strand_loci = find_overlapping_species_loci(proj, loci_for_species)
        if same_strand_loci:
            best, alternatives = choose_best_locus(proj, same_strand_loci)

            if best is not None:
                clocus.set_primary(proj.species, best.locus_id)

            if alternatives:
                clocus.set_alternatives(
                    proj.species,
                    [a.locus_id for a in alternatives]
                )
            continue

        # 2. any-strand overlaps: strand conflict candidates
        any_strand_loci = find_overlapping_species_loci_any_strand(proj, loci_for_species)
        if any_strand_loci:
            for locus in any_strand_loci:
                if locus.strand != proj.strand:
                    clocus.add_strand_conflict(proj.species, locus.locus_id)
            # if there was any overlap at all, do not classify as missing
            if proj.species in clocus.strand_conflicts:
                continue

        # 3. no overlaps at all: missing annotation candidate
        clocus.add_missing_projection(
            proj.species,
            f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}"
        )

    return clocus
