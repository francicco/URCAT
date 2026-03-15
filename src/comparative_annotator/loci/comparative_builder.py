from __future__ import annotations

from comparative_annotator.models.comparative import ComparativeLocus
from comparative_annotator.projection.matching import (
    match_projected_transcript_to_loci,
    classify_unmatched_projection,
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

        loci = species_loci.get(proj.species, [])

        matches = match_projected_transcript_to_loci(proj, loci)

        if matches:

            if matches[0].classification == "ortholog_candidate":

                clocus.add_ortholog(proj.species, matches[0].locus_id)

            else:

                for m in matches:
                    clocus.add_paralog(proj.species, m.locus_id)

        else:

            cls = classify_unmatched_projection(proj, loci)

            if cls == "missing_annotation_candidate":

                clocus.add_missing_projection(
                    proj.species,
                    f"{proj.seqid}:{proj.start}-{proj.end}",
                )

    return clocus
