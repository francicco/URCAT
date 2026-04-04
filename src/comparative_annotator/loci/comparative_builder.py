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
):
    """
    Build a ComparativeLocus object from projected transcripts.

    Logic:
    -------
    For each projected transcript:
        - find overlapping loci
        - classify:
            * primary (best overlap, same strand)
            * alternatives (other same-strand overlaps)
            * strand conflicts (overlap but wrong strand)
            * missing (no overlaps)

    Additionally:
        - if NO projected transcripts exist → record missing
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
        # This is CRITICAL — previously missing
        clocus.add_missing_projection(
            species="UNKNOWN",  # fallback, pipeline should override if needed
            projection_id=f"NO_PROJECTION:{seed_transcript}",
        )
        return clocus

    # --------------------------------------------------
    # Process each projected transcript
    # --------------------------------------------------
    for proj in projected_transcripts:
        target_species = proj.species

        loci = species_loci.get(target_species, [])
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
                target_species,
                f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}",
            )
            continue

        # --------------------------------------------------
        # Partition by strand
        # --------------------------------------------------
        same_strand = []
        opposite_strand = []

        for locus, frac in overlapping:
            if locus.strand == proj.strand:
                same_strand.append((locus, frac))
            else:
                opposite_strand.append((locus, frac))

        # --------------------------------------------------
        # Strand conflicts
        # --------------------------------------------------
        for locus, _ in opposite_strand:
            clocus.add_strand_conflict(target_species, locus.locus_id)

        # --------------------------------------------------
        # No same-strand → only conflicts → treat as missing
        # --------------------------------------------------
        if not same_strand:
            clocus.add_missing_projection(
                target_species,
                f"{proj.seqid}:{proj.start}-{proj.end}:{proj.strand}",
            )
            continue

        # --------------------------------------------------
        # Select primary (best overlap)
        # --------------------------------------------------
        same_strand.sort(key=lambda x: x[1], reverse=True)

        best_locus, best_score = same_strand[0]
        clocus.set_primary(target_species, best_locus.locus_id)

        # --------------------------------------------------
        # Set primary transcript if available
        # --------------------------------------------------
        if transcripts_by_species:
            txs = transcripts_by_species.get(target_species, {})
            if best_locus.locus_id in txs:
                clocus.set_primary_transcript(
                    target_species,
                    txs[best_locus.locus_id],
                )

        # --------------------------------------------------
        # Alternatives (remaining same-strand loci)
        # --------------------------------------------------
        alt_ids = [l.locus_id for l, _ in same_strand[1:]]

        if alt_ids:
            clocus.set_alternatives(target_species, alt_ids)

            if transcripts_by_species:
                txs = transcripts_by_species.get(target_species, {})
                alt_txs = [txs[lid] for lid in alt_ids if lid in txs]
                if alt_txs:
                    clocus.set_alternative_transcripts(target_species, alt_txs)

    return clocus
