from __future__ import annotations

from comparative_annotator.models.comparative import ComparativeLocus
from comparative_annotator.projection.adjudication import choose_best_locus
from comparative_annotator.projection.matching import (
    find_overlapping_species_loci,
    find_overlapping_species_loci_any_strand,
)
from comparative_annotator.projection.transcript_ranking import (
    choose_best_transcript_within_locus,
)


def build_comparative_locus_from_projection(
    seed_species,
    seed_transcript,
    projected_transcripts,
    species_loci,
    source_transcript,
    transcripts_by_species=None,
):
    if transcripts_by_species is None:
        transcripts_by_species = {}

    clocus = ComparativeLocus(
        locus_id=f"{seed_species}:{seed_transcript}",
        seed_species=seed_species,
        seed_transcript=seed_transcript,
    )

    for proj in projected_transcripts:
        loci_for_species = species_loci.get(proj.species, [])

        same_strand_loci = find_overlapping_species_loci(proj, loci_for_species)
        if same_strand_loci:
            best, alternatives = choose_best_locus(
                projected=proj,
                source_transcript=source_transcript,
                loci=same_strand_loci,
            )

            if best is not None:
                clocus.set_primary(proj.species, best.locus_id)

                tx_lookup = transcripts_by_species.get(proj.species, {})
                best_locus_obj = next((l for l in same_strand_loci if l.locus_id == best.locus_id), None)

                if best_locus_obj is not None and tx_lookup:
                    best_tx, alt_txs = choose_best_transcript_within_locus(
                        projected=proj,
                        source_transcript=source_transcript,
                        locus=best_locus_obj,
                        target_transcripts_by_id=tx_lookup,
                    )

                    if best_tx is not None:
                        clocus.set_primary_transcript(proj.species, best_tx.transcript_id)

                    if alt_txs:
                        clocus.set_alternative_transcripts(
                            proj.species,
                            [x.transcript_id for x in alt_txs],
                        )

            if alternatives:
                clocus.set_alternatives(
                    proj.species,
                    [a.locus_id for a in alternatives]
                )
            continue

        any_strand_loci = find_overlapping_species_loci_any_strand(proj, loci_for_species)
        if any_strand_loci:
            for locus in any_strand_loci:
                if locus.strand != proj.strand:
                    clocus.add_strand_conflict(proj.species, locus.locus_id)
            if proj.species in clocus.strand_conflicts:
                continue

        clocus.add_missing_projection(
            proj.species,
            f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}"
        )

    return clocus
