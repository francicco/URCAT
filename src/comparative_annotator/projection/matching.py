from __future__ import annotations

from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.match import ProjectionLocusMatch


def overlaps(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


def find_overlapping_species_loci(
    projected: ProjectedTranscript,
    species_loci: list[SpeciesLocus],
):
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
    overlap_start = max(projected.start, locus.start)
    overlap_end = min(projected.end, locus.end)

    overlap = max(0, overlap_end - overlap_start)

    proj_len = projected.end - projected.start

    if proj_len == 0:
        return 0.0

    return overlap / proj_len


def projected_transcript_id(projected: ProjectedTranscript) -> str:
    return (
        f"{projected.source_species}:{projected.source_transcript}"
        f"->{projected.species}:{projected.seqid}:{projected.start}-{projected.end}:{projected.strand}"
    )


def match_projected_transcript_to_loci(
    projected: ProjectedTranscript,
    species_loci: list[SpeciesLocus],
):
    overlapping = find_overlapping_species_loci(projected, species_loci)

    if not overlapping:
        return []

    classification = "ortholog_candidate" if len(overlapping) == 1 else "paralog_candidate"

    matches = []
    pid = projected_transcript_id(projected)

    for locus in overlapping:
        matches.append(
            ProjectionLocusMatch(
                projected_id=pid,
                species=projected.species,
                locus_id=locus.locus_id,
                overlap_fraction=locus_overlap_fraction(projected, locus),
                classification=classification,
            )
        )

    return sorted(matches, key=lambda x: x.overlap_fraction, reverse=True)


def classify_unmatched_projection(
    projected: ProjectedTranscript,
    species_loci: list[SpeciesLocus],
):
    overlapping = find_overlapping_species_loci(projected, species_loci)
    if overlapping:
        return None
    return "missing_annotation_candidate"


def nearest_species_locus(
    projected: ProjectedTranscript,
    species_loci: list[SpeciesLocus],
):
    """
    Return the nearest locus on the same seqid, regardless of overlap.
    """
    candidates = [l for l in species_loci if l.seqid == projected.seqid]
    if not candidates:
        return None, None

    best_locus = None
    best_distance = None

    for locus in candidates:
        if projected.end < locus.start:
            dist = locus.start - projected.end
        elif locus.end < projected.start:
            dist = projected.start - locus.end
        else:
            dist = 0

        if best_distance is None or dist < best_distance:
            best_distance = dist
            best_locus = locus

    return best_locus, best_distance


def find_overlapping_species_loci_any_strand(projected, species_loci):
    hits = []

    for locus in species_loci:
        if locus.seqid != projected.seqid:
            continue

        if overlaps(projected.start, projected.end, locus.start, locus.end):
            hits.append(locus)

    return hits
