from __future__ import annotations

from comparative_annotator.models.comparative import ComparativeLocus
from comparative_annotator.projection.adjudication import choose_best_locus
from comparative_annotator.projection.matching import find_overlapping_species_loci


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
        loci = find_overlapping_species_loci(
            proj,
            species_loci.get(proj.species, [])
        )

        if not loci:
            clocus.add_missing_projection(
                proj.species,
                f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}"
            )
            continue

        best, alternatives = choose_best_locus(proj, loci)

        if best is not None:
            clocus.set_primary(proj.species, best.locus_id)

        if alternatives:
            clocus.set_alternatives(
                proj.species,
                [a.locus_id for a in alternatives]
            )

    return clocus
