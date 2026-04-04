from __future__ import annotations

from comparative_annotator.models.comparative import ComparativeLocus


def _overlap_fraction(proj, locus) -> float:
    """
    Compute overlap fraction between projected transcript span and locus span.
    Uses 1-based inclusive coordinates.
    """
    p_start, p_end = proj.start, proj.end

    overlap_start = max(p_start, locus.start)
    overlap_end = min(p_end, locus.end)

    overlap = max(0, overlap_end - overlap_start + 1)
    proj_len = p_end - p_start + 1

    if proj_len <= 0:
        return 0.0

    return overlap / proj_len


def build_comparative_locus_from_projection(
    seed_species,
    seed_transcript,
    projected_transcripts,
    species_loci,
    source_transcript=None,
    transcripts_by_species=None,
    target_species=None,
):
    """
    Build a ComparativeLocus object from projected transcripts.

    Logic:
    -------
    For each projected transcript:
        - find overlapping loci
        - classify:
            * primary
            * alternatives
            * strand conflicts
            * missing

    Additionally:
        - if NO projected transcripts exist, record a missing projection for the
          requested target_species.
    """

    clocus = ComparativeLocus(
        locus_id=f"{seed_species}:{seed_transcript}",
        seed_species=seed_species,
        seed_transcript=seed_transcript,
    )

    # --------------------------------------------------
    # CASE 1: no projected transcripts at all
    # --------------------------------------------------
    if not projected_transcripts:
        if target_species is None:
            raise ValueError(
                "target_species must be provided when projected_transcripts is empty"
            )
        clocus.add_missing_projection(
            target_species,
            f"NO_PROJECTION:{seed_transcript}",
        )
        return clocus

    # --------------------------------------------------
    # Process each projected transcript
    # --------------------------------------------------
    for proj in projected_transcripts:
        proj_target_species = proj.species
        loci = species_loci.get(proj_target_species, [])
        overlapping = []

        for locus in loci:
            frac = _overlap_fraction(proj, locus)
            if frac > 0:
                overlapping.append((locus, frac))

        # --------------------------------------------------
        # CASE 2: no overlapping locus → missing
        # --------------------------------------------------
        if not overlapping:
            clocus.add_missing_projection(
                proj_target_species,
                f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}",
            )
            continue

        same_strand = []
        opposite_strand = []

        for locus, frac in overlapping:
            if locus.strand == proj.strand:
                same_strand.append((locus, frac))
            else:
                opposite_strand.append((locus, frac))

        for locus, _ in opposite_strand:
            clocus.add_strand_conflict(proj_target_species, locus.locus_id)

        # only wrong-strand overlaps -> still missing annotation
        if not same_strand:
            clocus.add_missing_projection(
                proj_target_species,
                f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}",
            )
            continue

        same_strand.sort(key=lambda x: x[1], reverse=True)

        best_locus, _ = same_strand[0]
        clocus.set_primary(proj_target_species, best_locus.locus_id)

        if transcripts_by_species:
            txs = transcripts_by_species.get(proj_target_species, {})
            best_tx = txs.get(best_locus.locus_id)
            if best_tx:
                if isinstance(best_tx, str):
                    clocus.set_primary_transcript(proj_target_species, best_tx)
                else:
                    tx_id = getattr(best_tx, "transcript_id", None)
                    if tx_id:
                        clocus.set_primary_transcript(proj_target_species, tx_id)

        alt_ids = [l.locus_id for l, _ in same_strand[1:]]
        if alt_ids:
            clocus.set_alternatives(proj_target_species, alt_ids)

            if transcripts_by_species:
                txs = transcripts_by_species.get(proj_target_species, {})
                alt_txs = []
                for lid in alt_ids:
                    tx_obj = txs.get(lid)
                    if tx_obj is None:
                        continue
                    if isinstance(tx_obj, str):
                        alt_txs.append(tx_obj)
                    else:
                        tx_id = getattr(tx_obj, "transcript_id", None)
                        if tx_id:
                            alt_txs.append(tx_id)
                if alt_txs:
                    clocus.set_alternative_transcripts(proj_target_species, alt_txs)

    return clocus
