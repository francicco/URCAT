from __future__ import annotations

from dataclasses import dataclass, field
from collections import defaultdict
from typing import Any

from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.transcript import CandidateTranscript


@dataclass(frozen=True)
class Interval:
    seqid: str
    start: int
    end: int
    strand: str | None = None

    def __post_init__(self) -> None:
        if self.end < self.start:
            raise ValueError(f"Invalid interval: {self.start}>{self.end}")

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Neighborhood:
    center_locus_id: str
    species: str
    seqid: str
    center_index: int
    left_neighbors: list[str] = field(default_factory=list)
    right_neighbors: list[str] = field(default_factory=list)
    window_neighbors: list[str] = field(default_factory=list)


@dataclass
class SpeciesIndex:
    species: str
    loci_by_seqid: dict[str, list[SpeciesLocus]]
    locus_id_to_locus: dict[str, SpeciesLocus]
    locus_id_to_position: dict[str, tuple[str, int]]
    neighborhood_cache: dict[str, Neighborhood] = field(default_factory=dict)


@dataclass
class AnchorMap:
    locus_to_anchor: dict[str, str] = field(default_factory=dict)

    def resolve(self, locus_id: str) -> str | None:
        return self.locus_to_anchor.get(locus_id)


@dataclass
class CandidateEdge:
    source_species: str
    source_locus_id: str
    target_species: str
    target_locus_id: str
    edge_origin: str
    projected_interval: Interval | None = None

    projection_features: dict[str, Any] = field(default_factory=dict)
    sequence_features: dict[str, Any] = field(default_factory=dict)
    architecture_features: dict[str, Any] = field(default_factory=dict)
    synteny_features: dict[str, Any] = field(default_factory=dict)
    flags: dict[str, Any] = field(default_factory=dict)

    edge_class: str | None = None
    edge_confidence: str | None = None
    accepted: bool | None = None

    def to_row(self) -> dict[str, Any]:
        row = {
            "source_species": self.source_species,
            "source_locus_id": self.source_locus_id,
            "target_species": self.target_species,
            "target_locus_id": self.target_locus_id,
            "edge_origin": self.edge_origin,
            "edge_class": self.edge_class,
            "edge_confidence": self.edge_confidence,
            "accepted": self.accepted,
        }
        row.update({f"proj_{k}": v for k, v in self.projection_features.items()})
        row.update({f"seq_{k}": v for k, v in self.sequence_features.items()})
        row.update({f"arch_{k}": v for k, v in self.architecture_features.items()})
        row.update({f"syn_{k}": v for k, v in self.synteny_features.items()})
        row.update({f"flag_{k}": v for k, v in self.flags.items()})
        return row


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return not (a_end < b_start or b_end < a_start)


def interval_overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if not interval_overlap(a_start, a_end, b_start, b_end):
        return 0
    return min(a_end, b_end) - max(a_start, b_start) + 1


def interval_distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if interval_overlap(a_start, a_end, b_start, b_end):
        return 0
    if a_end < b_start:
        return b_start - a_end
    return a_start - b_end


def jaccard_sets(a: set[str], b: set[str]) -> float | None:
    union = a | b
    if not union:
        return None
    return len(a & b) / len(union)


def _simple_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    if L == 0:
        return 0.0
    matches = sum(1 for i in range(L) if a[i] == b[i])
    return matches / L


class OrthologyIndexer:
    """
    Per-species locus ordering and neighborhood lookup.
    This should be created once per worker/job and reused.
    """

    def __init__(self, loci_by_species: dict[str, list[SpeciesLocus]]) -> None:
        self.species_indexes: dict[str, SpeciesIndex] = {}
        for species, loci in loci_by_species.items():
            self.species_indexes[species] = self._build_species_index(species, loci)

    @staticmethod
    def _build_species_index(species: str, loci: list[SpeciesLocus]) -> SpeciesIndex:
        loci_by_seqid: dict[str, list[SpeciesLocus]] = defaultdict(list)
        locus_id_to_locus: dict[str, SpeciesLocus] = {}
        locus_id_to_position: dict[str, tuple[str, int]] = {}

        for locus in loci:
            loci_by_seqid[locus.seqid].append(locus)
            locus_id_to_locus[locus.locus_id] = locus

        for seqid, seq_loci in loci_by_seqid.items():
            seq_loci.sort(key=lambda x: (x.start, x.end, x.locus_id))
            for i, locus in enumerate(seq_loci):
                locus_id_to_position[locus.locus_id] = (seqid, i)

        return SpeciesIndex(
            species=species,
            loci_by_seqid=dict(loci_by_seqid),
            locus_id_to_locus=locus_id_to_locus,
            locus_id_to_position=locus_id_to_position,
        )

    def get_species_index(self, species: str) -> SpeciesIndex:
        return self.species_indexes[species]

    def compute_neighborhood(
        self,
        species: str,
        locus_id: str,
        k_neighbors: int = 3,
        flank_bp: int = 100_000,
    ) -> Neighborhood:
        idx = self.get_species_index(species)
        if locus_id in idx.neighborhood_cache:
            return idx.neighborhood_cache[locus_id]

        locus = idx.locus_id_to_locus[locus_id]
        seqid, center_i = idx.locus_id_to_position[locus_id]
        seq_loci = idx.loci_by_seqid[seqid]

        left_i = max(0, center_i - k_neighbors)
        right_i = min(len(seq_loci), center_i + k_neighbors + 1)

        left_neighbors = [x.locus_id for x in seq_loci[left_i:center_i]]
        right_neighbors = [x.locus_id for x in seq_loci[center_i + 1:right_i]]

        win_start = max(1, locus.start - flank_bp)
        win_end = locus.end + flank_bp
        window_neighbors = [
            other.locus_id
            for other in seq_loci
            if other.locus_id != locus_id and interval_overlap(other.start, other.end, win_start, win_end)
        ]

        nb = Neighborhood(
            center_locus_id=locus_id,
            species=species,
            seqid=seqid,
            center_index=center_i,
            left_neighbors=left_neighbors,
            right_neighbors=right_neighbors,
            window_neighbors=window_neighbors,
        )
        idx.neighborhood_cache[locus_id] = nb
        return nb


class FeatureComputer:
    """
    Feature computation for orthology edges.
    """

    def __init__(
        self,
        loci_by_species: dict[str, list[SpeciesLocus]],
        transcripts_by_species: dict[str, dict[str, CandidateTranscript]] | None = None,
        anchor_map: AnchorMap | None = None,
        sequences_by_species: dict[str, dict[str, dict[str, str]]] | None = None,
    ) -> None:
        self.indexer = OrthologyIndexer(loci_by_species)
        self.transcripts_by_species = transcripts_by_species or {}
        self.anchor_map = anchor_map or AnchorMap()
        self.sequences_by_species = sequences_by_species or {}
        self.diamond_cache = {}

    def build_edge(
        self,
        source_species: str,
        source_locus_id: str,
        target_species: str,
        target_locus_id: str,
        edge_origin: str,
        projected_interval: Interval | None = None,
    ) -> CandidateEdge:
        source_locus = self.indexer.get_species_index(source_species).locus_id_to_locus[source_locus_id]
        target_locus = self.indexer.get_species_index(target_species).locus_id_to_locus[target_locus_id]

        edge = CandidateEdge(
            source_species=source_species,
            source_locus_id=source_locus_id,
            target_species=target_species,
            target_locus_id=target_locus_id,
            edge_origin=edge_origin,
            projected_interval=projected_interval,
        )

        edge.projection_features = self.compute_projection_features(
            source_locus=source_locus,
            target_locus=target_locus,
            projected_interval=projected_interval,
        )
        edge.sequence_features = self.compute_sequence_features(
            source_species=source_species,
            source_locus=source_locus,
            target_species=target_species,
            target_locus=target_locus,
        )
        edge.architecture_features = self.compute_architecture_features(
            source_species=source_species,
            source_locus=source_locus,
            target_species=target_species,
            target_locus=target_locus,
        )
        edge.synteny_features = self.compute_microsynteny_features(
            source_species=source_species,
            source_locus=source_locus,
            target_species=target_species,
            target_locus=target_locus,
            projected_interval=projected_interval,
        )
        edge.flags = self.compute_edge_flags(edge)
        return edge

    def compute_projection_features(
        self,
        source_locus: SpeciesLocus,
        target_locus: SpeciesLocus,
        projected_interval: Interval | None,
    ) -> dict[str, Any]:
        proj_cov = None
        anchor_distance = None

        if projected_interval is not None and projected_interval.seqid == target_locus.seqid:
            overlap = interval_overlap_len(
                target_locus.start,
                target_locus.end,
                projected_interval.start,
                projected_interval.end,
            )
            denom = max(1, min(target_locus.end - target_locus.start + 1, projected_interval.length))
            proj_cov = overlap / denom
            anchor_distance = interval_distance(
                target_locus.start,
                target_locus.end,
                projected_interval.start,
                projected_interval.end,
            )

        return {
            "cov": proj_cov,
            "same_strand": (
                projected_interval is None
                or projected_interval.strand is None
                or projected_interval.strand == target_locus.strand
            ),
            "anchor_distance": anchor_distance,
            "projected_seqid_match": (
                projected_interval is None or projected_interval.seqid == target_locus.seqid
            ),
        }

    def compute_sequence_features(
        self,
        source_species: str,
        source_locus: SpeciesLocus,
        target_species: str,
        target_locus: SpeciesLocus,
    ) -> dict[str, Any]:
        src_txs = self._get_locus_transcripts(source_species, source_locus)
        tgt_txs = self._get_locus_transcripts(target_species, target_locus)

        src_aa = self.sequences_by_species.get(source_species, {}).get("aa", {})
        hits = self.diamond_cache.get((source_species, target_species), {})

        best_score = -1.0
        best_pair = (None, None)
        best_hit = None

        for s in src_txs:
            for t in tgt_txs:
                key = (s.transcript_id, t.transcript_id)
                hit = hits.get(key)
                if hit is None:
                    continue

                score = hit["bitscore"]
                if score > best_score:
                    best_score = score
                    best_pair = (s.transcript_id, t.transcript_id)
                    best_hit = hit

        if best_hit is None:
            return {
                "prot_cov": None,
                "prot_id": None,
                "bitscore": None,
                "best_hit_margin": None,
                "cds_intact": None,
                "cds_complete": None,
                "best_source_tx": None,
                "best_target_tx": None,
            }

        src_len = len(src_aa.get(best_pair[0], ""))

        return {
            "prot_cov": best_hit["aln_len"] / max(1, src_len),
            "prot_id": best_hit["pid"],
            "bitscore": best_hit["bitscore"],
            "best_hit_margin": None,
            "cds_intact": True,
            "cds_complete": True,
            "best_source_tx": best_pair[0],
            "best_target_tx": best_pair[1],
        }

    def compute_architecture_features(
        self,
        source_species: str,
        source_locus: SpeciesLocus,
        target_species: str,
        target_locus: SpeciesLocus,
    ) -> dict[str, Any]:
        src_txs = self._get_locus_transcripts(source_species, source_locus)
        tgt_txs = self._get_locus_transcripts(target_species, target_locus)

        src_introns = self._collect_intron_signatures(src_txs)
        tgt_introns = self._collect_intron_signatures(tgt_txs)

        shared_introns = len(src_introns & tgt_introns)
        intron_recall = None if not src_introns else shared_introns / len(src_introns)
        intron_precision = None if not tgt_introns else shared_introns / len(tgt_introns)

        src_best_exons = max((tx.exon_count for tx in src_txs), default=0)
        tgt_best_exons = max((tx.exon_count for tx in tgt_txs), default=0)

        return {
            "shared_introns": shared_introns,
            "intron_recall": intron_recall,
            "intron_precision": intron_precision,
            "exon_count_compatible": abs(src_best_exons - tgt_best_exons) <= 1,
            "utr_support": None,
        }

    def compute_microsynteny_features(
        self,
        source_species: str,
        source_locus: SpeciesLocus,
        target_species: str,
        target_locus: SpeciesLocus,
        projected_interval: Interval | None,
        k_neighbors: int = 3,
        flank_bp: int = 100_000,
        anchor_tolerance_bp: int = 50_000,
    ) -> dict[str, Any]:
        source_nb = self.indexer.compute_neighborhood(source_species, source_locus.locus_id, k_neighbors, flank_bp)
        target_nb = self.indexer.compute_neighborhood(target_species, target_locus.locus_id, k_neighbors, flank_bp)

        anchor_ok = self._is_candidate_in_expected_anchor(
            candidate_target=target_locus,
            projected_interval=projected_interval,
            anchor_tolerance_bp=anchor_tolerance_bp,
        )
        anchor_distance = self._distance_to_projected_anchor(
            candidate_target=target_locus,
            projected_interval=projected_interval,
        )

        source_left_anchor, source_right_anchor = self._nearest_left_right_anchor_ids(source_nb)
        target_left_anchor, target_right_anchor = self._nearest_left_right_anchor_ids(target_nb)

        left_neighbor_match = int(
            source_left_anchor is not None
            and target_left_anchor is not None
            and source_left_anchor == target_left_anchor
        )
        right_neighbor_match = int(
            source_right_anchor is not None
            and target_right_anchor is not None
            and source_right_anchor == target_right_anchor
        )

        source_anchor_set = self._neighborhood_anchor_set(source_nb)
        target_anchor_set = self._neighborhood_anchor_set(target_nb)
        neighbor_jaccard = jaccard_sets(source_anchor_set, target_anchor_set)
        order_score = self._compute_local_order_score(source_nb, target_nb)
        orientation_score = self._compute_local_orientation_score(source_nb, target_nb)

        neighbor_hits = left_neighbor_match + right_neighbor_match

        synteny_class = self._classify_synteny(
            anchor_ok=anchor_ok,
            neighbor_hits=neighbor_hits,
            neighbor_jaccard=neighbor_jaccard,
            order_score=order_score,
        )

        return {
            "anchor_ok": anchor_ok,
            "anchor_distance": anchor_distance,
            "left_neighbor_match": left_neighbor_match,
            "right_neighbor_match": right_neighbor_match,
            "neighbor_hits": neighbor_hits,
            "neighbor_jaccard": neighbor_jaccard,
            "order_score": order_score,
            "orientation_score": orientation_score,
            "synteny_class": synteny_class,
        }

    def compute_edge_flags(self, edge: CandidateEdge) -> dict[str, Any]:
        return {
            "strand_conflict": edge.projection_features.get("same_strand") is False,
            "has_competitor": False,
        }

    def _get_locus_transcripts(self, species: str, locus: SpeciesLocus) -> list[CandidateTranscript]:
        tx_by_id = self.transcripts_by_species.get(species, {})
        return [tx_by_id[txid] for txid in locus.transcripts if txid in tx_by_id]

    @staticmethod
    def _collect_intron_signatures(transcripts: list[CandidateTranscript]) -> set[tuple[int, int]]:
        out: set[tuple[int, int]] = set()
        for tx in transcripts:
            out.update(tx.intron_chain)
        return out

    @staticmethod
    def _is_candidate_in_expected_anchor(
        candidate_target: SpeciesLocus,
        projected_interval: Interval | None,
        anchor_tolerance_bp: int,
    ) -> bool:
        if projected_interval is None:
            return False
        if candidate_target.seqid != projected_interval.seqid:
            return False
        return interval_overlap(
            candidate_target.start,
            candidate_target.end,
            max(1, projected_interval.start - anchor_tolerance_bp),
            projected_interval.end + anchor_tolerance_bp,
        )

    @staticmethod
    def _distance_to_projected_anchor(
        candidate_target: SpeciesLocus,
        projected_interval: Interval | None,
    ) -> int | None:
        if projected_interval is None:
            return None
        if candidate_target.seqid != projected_interval.seqid:
            return None
        return interval_distance(
            candidate_target.start,
            candidate_target.end,
            projected_interval.start,
            projected_interval.end,
        )

    def _neighborhood_anchor_set(self, nb: Neighborhood) -> set[str]:
        anchors = set()
        for locus_id in nb.window_neighbors:
            aid = self.anchor_map.resolve(locus_id)
            if aid is not None:
                anchors.add(aid)
        return anchors

    def _nearest_left_right_anchor_ids(self, nb: Neighborhood) -> tuple[str | None, str | None]:
        left_anchor = None
        right_anchor = None

        for locus_id in reversed(nb.left_neighbors):
            aid = self.anchor_map.resolve(locus_id)
            if aid is not None:
                left_anchor = aid
                break

        for locus_id in nb.right_neighbors:
            aid = self.anchor_map.resolve(locus_id)
            if aid is not None:
                right_anchor = aid
                break

        return left_anchor, right_anchor

    def _compute_local_order_score(self, source_nb: Neighborhood, target_nb: Neighborhood) -> float | None:
        source_order = [self.anchor_map.resolve(x) for x in source_nb.window_neighbors]
        target_order = [self.anchor_map.resolve(x) for x in target_nb.window_neighbors]
        source_order = [x for x in source_order if x is not None]
        target_order = [x for x in target_order if x is not None]

        shared = [x for x in source_order if x in set(target_order)]
        if len(shared) < 2:
            return None

        src_rank = {a: i for i, a in enumerate(source_order)}
        tgt_rank = {a: i for i, a in enumerate(target_order)}

        concordant = 0
        total = 0
        for i in range(len(shared)):
            for j in range(i + 1, len(shared)):
                a = shared[i]
                b = shared[j]
                total += 1
                if (src_rank[a] < src_rank[b]) == (tgt_rank[a] < tgt_rank[b]):
                    concordant += 1

        return None if total == 0 else concordant / total

    def _compute_local_orientation_score(self, source_nb: Neighborhood, target_nb: Neighborhood) -> float | None:
        source_idx = self.indexer.get_species_index(source_nb.species)
        target_idx = self.indexer.get_species_index(target_nb.species)

        src_map: dict[str, str] = {}
        tgt_map: dict[str, str] = {}

        for locus_id in source_nb.window_neighbors:
            aid = self.anchor_map.resolve(locus_id)
            if aid is not None:
                src_map[aid] = source_idx.locus_id_to_locus[locus_id].strand

        for locus_id in target_nb.window_neighbors:
            aid = self.anchor_map.resolve(locus_id)
            if aid is not None:
                tgt_map[aid] = target_idx.locus_id_to_locus[locus_id].strand

        shared = set(src_map) & set(tgt_map)
        if not shared:
            return None

        matches = sum(1 for aid in shared if src_map[aid] == tgt_map[aid])
        return matches / len(shared)

    @staticmethod
    def _classify_synteny(
        anchor_ok: bool,
        neighbor_hits: int,
        neighbor_jaccard: float | None,
        order_score: float | None,
    ) -> str:
        j = neighbor_jaccard or 0.0
        o = order_score or 0.0

        if anchor_ok and (neighbor_hits == 2 or j >= 0.5):
            return "strong"
        if anchor_ok and (neighbor_hits >= 1 or j >= 0.25):
            return "moderate"
        if anchor_ok:
            return "weak"
        if (neighbor_hits >= 1 and j >= 0.25) or o >= 0.5:
            return "partial_nonanchor"
        return "none"


class RuleBasedEdgeClassifier:
    def classify_edges(self, edges: list[CandidateEdge]) -> list[CandidateEdge]:
        grouped: dict[tuple[str, str], list[CandidateEdge]] = defaultdict(list)
        for edge in edges:
            grouped[(edge.source_species, edge.source_locus_id)].append(edge)

        out = []
        for _, comp_edges in grouped.items():
            for edge in comp_edges:
                others = [x for x in comp_edges if x.target_locus_id != edge.target_locus_id]
                edge.edge_class, edge.edge_confidence, edge.accepted = self.classify_edge(edge, others)
                out.append(edge)
        return out

    def classify_edge(
        self,
        edge: CandidateEdge,
        competing_edges: list[CandidateEdge],
    ) -> tuple[str, str, bool]:
        p = edge.projection_features
        s = edge.sequence_features
        a = edge.architecture_features
        y = edge.synteny_features

        projection_pass = self._passes_projection(p)
        sequence_pass = self._passes_sequence(s)
        architecture_pass = self._passes_architecture(a)
        synteny_pass = y.get("synteny_class") in {"strong", "moderate", "weak"}
        moderate_synteny = y.get("synteny_class") in {"strong", "moderate"}

        better_syntenic_competitor = self._has_better_syntenic_competitor(edge, competing_edges)

        if not projection_pass and not sequence_pass and y.get("synteny_class") == "none":
            return "N", "low", False

        if better_syntenic_competitor and sequence_pass and not moderate_synteny:
            return "P", "medium", False

        if projection_pass and moderate_synteny and architecture_pass and s.get("cds_intact") is False:
            return "E", "medium", True

        if self._is_fragmentation_pattern(edge, competing_edges):
            return "D", "medium", True

        if projection_pass and sequence_pass and architecture_pass and moderate_synteny:
            return "A", "high", True

        if projection_pass and architecture_pass and moderate_synteny and not bool(s.get("cds_complete")):
            return "B", "medium", True

        if projection_pass and sequence_pass and synteny_pass:
            return "C", "medium", True

        if sequence_pass and not synteny_pass:
            return "P", "low", False

        return "X", "low", False

    @staticmethod
    def _passes_projection(p: dict[str, Any]) -> bool:
        cov = p.get("cov")
        same_strand = p.get("same_strand")
        return (cov is not None and cov >= 0.3) and (same_strand is not False)

    @staticmethod
    def _passes_sequence(s: dict[str, Any]) -> bool:
        prot_cov = s.get("prot_cov")
        prot_id = s.get("prot_id")
        if prot_cov is None or prot_id is None:
            return False
        return prot_cov >= 0.4 and prot_id >= 0.4

    @staticmethod
    def _passes_architecture(a: dict[str, Any]) -> bool:
        recall = a.get("intron_recall")
        precision = a.get("intron_precision")
        exon_ok = a.get("exon_count_compatible")

        if recall is None and precision is None:
            return bool(exon_ok)

        return (
            (recall is not None and recall >= 0.3)
            or (precision is not None and precision >= 0.3)
            or bool(exon_ok)
        )

    @staticmethod
    def _synteny_rank(edge: CandidateEdge) -> int:
        cls = edge.synteny_features.get("synteny_class")
        return {
            "strong": 4,
            "moderate": 3,
            "weak": 2,
            "partial_nonanchor": 1,
            "none": 0,
        }.get(cls, 0)

    def _has_better_syntenic_competitor(
        self,
        edge: CandidateEdge,
        competing_edges: list[CandidateEdge],
    ) -> bool:
        my_rank = self._synteny_rank(edge)
        my_seq = edge.sequence_features.get("prot_cov") or 0.0

        for other in competing_edges:
            other_rank = self._synteny_rank(other)
            other_seq = other.sequence_features.get("prot_cov") or 0.0
            if other_rank > my_rank and other_seq >= my_seq:
                return True
        return False

    @staticmethod
    def _is_fragmentation_pattern(
        edge: CandidateEdge,
        competing_edges: list[CandidateEdge],
    ) -> bool:
        if edge.synteny_features.get("synteny_class") not in {"strong", "moderate"}:
            return False

        my_cov = edge.projection_features.get("cov") or 0.0
        if my_cov > 0.8:
            return False

        n_local = 0
        for other in competing_edges:
            if other.target_species != edge.target_species:
                continue
            if other.synteny_features.get("synteny_class") in {"strong", "moderate"}:
                other_cov = other.projection_features.get("cov") or 0.0
                if 0.1 <= other_cov <= 0.8:
                    n_local += 1

        return n_local >= 1


class SimpleGraph:
    def __init__(self) -> None:
        self.adj: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)

    def add_node(self, node: tuple[str, str]) -> None:
        _ = self.adj[node]

    def add_edge(self, u: tuple[str, str], v: tuple[str, str]) -> None:
        self.adj[u].add(v)
        self.adj[v].add(u)

    def connected_components(self) -> list[set[tuple[str, str]]]:
        seen = set()
        components = []

        for node in self.adj:
            if node in seen:
                continue

            stack = [node]
            comp = set()
            seen.add(node)

            while stack:
                cur = stack.pop()
                comp.add(cur)
                for nxt in self.adj[cur]:
                    if nxt not in seen:
                        seen.add(nxt)
                        stack.append(nxt)

            components.append(comp)

        return components


class OrthogroupBuilder:
    def build_graph(self, edges: list[CandidateEdge], accepted_classes: set[str]) -> SimpleGraph:
        g = SimpleGraph()
        for edge in edges:
            if edge.accepted and edge.edge_class in accepted_classes:
                u = (edge.source_species, edge.source_locus_id)
                v = (edge.target_species, edge.target_locus_id)
                g.add_node(u)
                g.add_node(v)
                g.add_edge(u, v)
        return g

    def compute_orthogroups(self, edges: list[CandidateEdge]) -> dict[str, list[set[tuple[str, str]]]]:
        strict = self.build_graph(edges, {"A"}).connected_components()
        standard = self.build_graph(edges, {"A", "B", "C", "E"}).connected_components()
        exploratory = self.build_graph(edges, {"A", "B", "C", "D", "E", "X"}).connected_components()
        return {
            "strict": strict,
            "standard": standard,
            "exploratory": exploratory,
        }


def build_and_classify_edges(
    loci_by_species: dict[str, list[SpeciesLocus]],
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
    candidate_pairs: list[tuple[str, str, str, str, str, Interval | None]],
    anchor_map: AnchorMap | None = None,
    sequences_by_species: dict[str, dict[str, dict[str, str]]] | None = None,
    diamond_cache: dict | None = None,
) -> tuple[list[CandidateEdge], dict[str, list[set[tuple[str, str]]]]]:
    fc = FeatureComputer(
        loci_by_species=loci_by_species,
        transcripts_by_species=transcripts_by_species,
        anchor_map=anchor_map,
        sequences_by_species=sequences_by_species,
    )
    fc.diamond_cache = diamond_cache or {}
    
    edges = []
    for (
        source_species,
        source_locus_id,
        target_species,
        target_locus_id,
        edge_origin,
        projected_interval,
    ) in candidate_pairs:
        edge = fc.build_edge(
            source_species=source_species,
            source_locus_id=source_locus_id,
            target_species=target_species,
            target_locus_id=target_locus_id,
            edge_origin=edge_origin,
            projected_interval=projected_interval,
        )
        edges.append(edge)

    classifier = RuleBasedEdgeClassifier()
    classified = classifier.classify_edges(edges)

    builder = OrthogroupBuilder()
    orthogroups = builder.compute_orthogroups(classified)
    return classified, orthogroups
