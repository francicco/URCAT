from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass


@dataclass
class MissingLocusConsensus:
    species: str
    seqid: str
    strand: str

    support_count: int
    total_chain_score: float
    mean_exon_recovery: float

    source_transcripts: list[str]
    exons: list[tuple[int, int]]


def _projected_span(projected):
    return projected.start, projected.end


def _overlaps_or_is_near(a_start, a_end, b_start, b_end, max_gap=0):
    if a_end < b_start:
        return (b_start - a_end) <= max_gap
    if b_end < a_start:
        return (a_start - b_end) <= max_gap
    return True


def cluster_projected_transcripts(projected_transcripts, max_gap=0):
    """
    Cluster projected transcripts by species, seqid, and overlapping/nearby span.
    Strand is ignored at this step so opposite-strand hypotheses can compete
    for the same missing locus.
    """
    if not projected_transcripts:
        return []

    items = sorted(
        projected_transcripts,
        key=lambda p: (p.species, p.seqid, p.start, p.end),
    )

    clusters = []
    current = [items[0]]
    cur_species = items[0].species
    cur_seqid = items[0].seqid
    cur_start = items[0].start
    cur_end = items[0].end

    for pt in items[1:]:
        same_target = (pt.species == cur_species and pt.seqid == cur_seqid)
        if same_target and _overlaps_or_is_near(cur_start, cur_end, pt.start, pt.end, max_gap=max_gap):
            current.append(pt)
            cur_start = min(cur_start, pt.start)
            cur_end = max(cur_end, pt.end)
        else:
            clusters.append(current)
            current = [pt]
            cur_species = pt.species
            cur_seqid = pt.seqid
            cur_start = pt.start
            cur_end = pt.end

    if current:
        clusters.append(current)

    return clusters


def choose_missing_locus_strand(projected_transcripts):
    """
    Choose the best-supported strand hypothesis among projected transcripts
    that belong to the same missing-locus cluster.
    """
    if not projected_transcripts:
        return None, {}

    by_strand = defaultdict(list)
    for pt in projected_transcripts:
        by_strand[pt.strand].append(pt)

    strand_scores = {}

    for strand, pts in by_strand.items():
        support_count = len(pts)
        total_chain_score = sum((pt.chain_score or 0.0) for pt in pts)

        exon_recoveries = []
        for pt in pts:
            if pt.exon_count > 0 and pt.coverage is not None:
                exon_recoveries.append(pt.coverage / pt.exon_count)
            else:
                exon_recoveries.append(0.0)

        mean_exon_recovery = sum(exon_recoveries) / len(exon_recoveries) if exon_recoveries else 0.0

        score = (
            0.5 * support_count +
            0.3 * total_chain_score +
            0.2 * mean_exon_recovery
        )

        strand_scores[strand] = {
            "support_count": support_count,
            "total_chain_score": total_chain_score,
            "mean_exon_recovery": mean_exon_recovery,
            "score": score,
            "projected_transcripts": pts,
        }

    best_strand = max(strand_scores, key=lambda s: strand_scores[s]["score"])
    return best_strand, strand_scores


def _median(values):
    values = sorted(values)
    n = len(values)
    if n == 0:
        return None
    mid = n // 2
    if n % 2 == 1:
        return values[mid]
    return int(round((values[mid - 1] + values[mid]) / 2))


def build_consensus_missing_transcript(projected_transcripts, chosen_strand):
    """
    Build a simple consensus transcript from projected transcripts on the chosen strand.

    Version 1:
    - group transcripts by exon count
    - choose the most common exon count
    - compute median exon boundaries by exon index
    """
    pts = [pt for pt in projected_transcripts if pt.strand == chosen_strand]
    if not pts:
        return None

    species = pts[0].species
    seqid = pts[0].seqid

    # choose the dominant exon count
    exon_count_groups = defaultdict(list)
    for pt in pts:
        exon_count_groups[pt.exon_count].append(pt)

    best_exon_count = max(
        exon_count_groups,
        key=lambda n: (
            len(exon_count_groups[n]),
            sum((pt.chain_score or 0.0) for pt in exon_count_groups[n]),
        ),
    )

    supporting_pts = exon_count_groups[best_exon_count]

    consensus_exons = []
    for exon_idx in range(best_exon_count):
        starts = [pt.exons[exon_idx][0] for pt in supporting_pts]
        ends = [pt.exons[exon_idx][1] for pt in supporting_pts]

        consensus_exons.append((_median(starts), _median(ends)))

    support_count = len(supporting_pts)
    total_chain_score = sum((pt.chain_score or 0.0) for pt in supporting_pts)

    exon_recoveries = []
    for pt in supporting_pts:
        if pt.exon_count > 0 and pt.coverage is not None:
            exon_recoveries.append(pt.coverage / pt.exon_count)
        else:
            exon_recoveries.append(0.0)

    mean_exon_recovery = sum(exon_recoveries) / len(exon_recoveries) if exon_recoveries else 0.0

    return MissingLocusConsensus(
        species=species,
        seqid=seqid,
        strand=chosen_strand,
        support_count=support_count,
        total_chain_score=total_chain_score,
        mean_exon_recovery=mean_exon_recovery,
        source_transcripts=[pt.source_transcript for pt in supporting_pts],
        exons=consensus_exons,
    )
