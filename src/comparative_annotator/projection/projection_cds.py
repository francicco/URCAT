from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG"}


@dataclass
class ProjectedCDSBlock:
    source_start: int
    source_end: int
    target_seqid: str
    target_start: int
    target_end: int
    target_strand: str | None


@dataclass
class ProjectedCDSModel:
    source_species: str
    source_transcript_id: str
    target_species: str
    target_gene_id: str | None
    target_transcript_id: str | None
    seqid: str | None
    strand: str | None
    projected_exon_intervals: list[tuple[int, int]] = field(default_factory=list)
    projected_cds_intervals: list[tuple[int, int]] = field(default_factory=list)
    all_projected_cds_blocks: list[ProjectedCDSBlock] = field(default_factory=list)
    projected_seqids: list[str] = field(default_factory=list)
    source_cds_length_nt: int | None = None
    target_cds_sequence: str | None = None
    source_cds_sequence: str | None = None


@dataclass
class CodingIntegrityReport:
    cds_length_nt: int
    cds_length_aa: int
    cds_recovery: float | None
    has_start_codon: bool | None
    has_stop_codon: bool | None
    internal_stop_count: int | None
    length_mod_3: int | None
    is_in_frame: bool | None
    phase_series: list[int] = field(default_factory=list)
    source_cds_length_nt: int | None = None
    target_cds_interval_count: int = 0
    classification_hint: str = "unknown"


@dataclass
class AssemblyFragmentationReport:
    is_fragmented_across_seqids: bool
    n_seqids: int
    seqids: list[str]
    classification_hint: str


@dataclass
class ProjectedLocusAssessment:
    coding: CodingIntegrityReport
    fragmentation: AssemblyFragmentationReport
    final_class: str
    final_reason: str


@dataclass
class TranscriptProjectionLike:
    """
    Minimal structural interface expected from reconstruct_projected_transcripts().
    """

    source_species: str
    source_transcript: str
    target_species: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]]


# ---------------------------------------------------------------------------
# Basic sequence helpers
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def normalize_seq(seq: str | None) -> str | None:
    if seq is None:
        return None
    return "".join(seq.split()).upper()


def interval_length(start: int, end: int) -> int:
    return abs(end - start) + 1


def normalize_interval(start: int, end: int) -> tuple[int, int]:
    return (start, end) if start <= end else (end, start)


def sum_interval_lengths(intervals: list[tuple[int, int]]) -> int:
    return sum(interval_length(s, e) for s, e in intervals)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []

    ordered = sorted(normalize_interval(s, e) for s, e in intervals)
    merged: list[tuple[int, int]] = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def cds_intervals_in_transcript_order(
    intervals: list[tuple[int, int]],
    strand: str | None,
) -> list[tuple[int, int]]:
    ordered = sorted(normalize_interval(s, e) for s, e in intervals)
    if strand == "-":
        return list(reversed(ordered))
    return ordered


def extract_subsequence(
    genome_by_seqid: dict[str, str],
    seqid: str,
    start: int,
    end: int,
) -> str:
    start, end = normalize_interval(start, end)
    seq = genome_by_seqid[seqid]
    return seq[start - 1:end]


def extract_spliced_sequence(
    genome_by_seqid: dict[str, str],
    seqid: str,
    intervals: list[tuple[int, int]],
    strand: str | None,
) -> str:
    ordered = cds_intervals_in_transcript_order(intervals, strand)
    seq = "".join(extract_subsequence(genome_by_seqid, seqid, s, e) for s, e in ordered)
    if strand == "-":
        seq = reverse_complement(seq)
    return normalize_seq(seq) or ""


# ---------------------------------------------------------------------------
# Translation helpers
# ---------------------------------------------------------------------------

def translate_cds(seq: str) -> str:
    seq = normalize_seq(seq) or ""
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    aa: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        aa.append(codon_table.get(seq[i:i + 3], "X"))
    return "".join(aa)


def count_internal_stops(aa: str) -> int:
    if not aa:
        return 0
    body = aa[:-1] if aa.endswith("*") else aa
    return body.count("*")


def cds_phase_series(intervals: list[tuple[int, int]], strand: str | None) -> list[int]:
    ordered = cds_intervals_in_transcript_order(intervals, strand)
    phases_tx_order: list[int] = []
    consumed = 0
    for start, end in ordered:
        phase = 0 if consumed == 0 else (3 - (consumed % 3)) % 3
        phases_tx_order.append(phase)
        consumed += interval_length(start, end)
    if strand == "-":
        return list(reversed(phases_tx_order))
    return phases_tx_order


# ---------------------------------------------------------------------------
# Transcript attribute readers
# ---------------------------------------------------------------------------

def _coerce_interval_list(value: Any) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    if not value:
        return out
    for item in value:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            out.append(normalize_interval(int(item[0]), int(item[1])))
    return sorted(out)


def get_transcript_cds_intervals(tx: Any) -> list[tuple[int, int]]:
    for attr in ("cds_parts", "cds_intervals", "cds", "CDS", "cds_blocks"):
        value = getattr(tx, attr, None)
        out = _coerce_interval_list(value)
        if out:
            return out
    return []


def get_transcript_exon_intervals(tx: Any) -> list[tuple[int, int]]:
    for attr in ("exons", "exon_intervals", "exon_blocks"):
        value = getattr(tx, attr, None)
        out = _coerce_interval_list(value)
        if out:
            return out
    return []


# ---------------------------------------------------------------------------
# Projection helpers
# ---------------------------------------------------------------------------

def project_intervals_with_hal(
    hal: Any,
    source_species: str,
    target_species: str,
    source_seqid: str,
    source_strand: str,
    intervals: list[tuple[int, int]],
    source_transcript_id: str | None = None,
) -> list[ProjectedCDSBlock]:
    out: list[ProjectedCDSBlock] = []
    for source_start, source_end in intervals:
        hits = hal.project_interval(
            source_species=source_species,
            target_species=target_species,
            seqid=source_seqid,
            start=source_start,
            end=source_end,
            strand=source_strand,
            source_transcript=source_transcript_id,
        )
        for h in hits:
            out.append(
                ProjectedCDSBlock(
                    source_start=int(source_start),
                    source_end=int(source_end),
                    target_seqid=h.seqid,
                    target_start=int(h.start),
                    target_end=int(h.end),
                    target_strand=h.strand if h.strand is not None else source_strand,
                )
            )
    return out


def _dominant_seqid_from_blocks(blocks: list[ProjectedCDSBlock]) -> str | None:
    if not blocks:
        return None
    support: dict[str, int] = {}
    for b in blocks:
        support[b.target_seqid] = support.get(b.target_seqid, 0) + interval_length(b.target_start, b.target_end)
    return sorted(support.items(), key=lambda x: (-x[1], x[0]))[0][0]


def _dominant_strand_from_blocks(blocks: list[ProjectedCDSBlock], default: str | None = None) -> str | None:
    if not blocks:
        return default
    support: dict[str, int] = {}
    for b in blocks:
        strand = b.target_strand or default or "."
        support[strand] = support.get(strand, 0) + 1
    return sorted(support.items(), key=lambda x: (-x[1], x[0]))[0][0]


def restrict_blocks_to_projected_transcript(
    projected_cds_blocks: list[ProjectedCDSBlock],
    projected_tx: TranscriptProjectionLike | None,
) -> tuple[list[ProjectedCDSBlock], list[str]]:
    """
    Keep only CDS blocks that are compatible with the accepted projected transcript.

    Rules:
    - if projected_tx is absent, fall back to dominant-seqid blocks
    - otherwise require same seqid and strand as projected_tx
    - then keep only CDS blocks that overlap at least one projected exon
    """
    if not projected_cds_blocks:
        return [], []

    all_seqids = sorted({b.target_seqid for b in projected_cds_blocks})

    if projected_tx is None:
        dominant_seqid = _dominant_seqid_from_blocks(projected_cds_blocks)
        kept = [b for b in projected_cds_blocks if b.target_seqid == dominant_seqid]
        return kept, all_seqids

    exon_intervals = [normalize_interval(s, e) for s, e in projected_tx.exons]
    tx_seqid = projected_tx.seqid
    tx_strand = projected_tx.strand

    kept: list[ProjectedCDSBlock] = []
    for b in projected_cds_blocks:
        if b.target_seqid != tx_seqid:
            continue
        if tx_strand is not None and b.target_strand is not None and b.target_strand != tx_strand:
            continue
        bs, be = normalize_interval(b.target_start, b.target_end)
        overlaps_exon = any(not (be < es or bs > ee) for es, ee in exon_intervals)
        if overlaps_exon:
            kept.append(b)
    return kept, all_seqids


def collapse_blocks_to_intervals(blocks: list[ProjectedCDSBlock]) -> list[tuple[int, int]]:
    return merge_intervals([(b.target_start, b.target_end) for b in blocks])


# ---------------------------------------------------------------------------
# Main model builder
# ---------------------------------------------------------------------------

def build_projected_cds_model(
    source_tx: Any,
    projected_tx: TranscriptProjectionLike | None,
    hal: Any,
    target_species: str,
    target_gene_id: str | None = None,
    target_transcript_id: str | None = None,
    target_genome_by_seqid: dict[str, str] | None = None,
    source_genome_by_seqid: dict[str, str] | None = None,
) -> ProjectedCDSModel:
    source_cds = get_transcript_cds_intervals(source_tx)
    source_cds_len = sum_interval_lengths(source_cds) if source_cds else None

    projected_blocks = project_intervals_with_hal(
        hal=hal,
        source_species=source_tx.species,
        target_species=target_species,
        source_seqid=source_tx.seqid,
        source_strand=source_tx.strand,
        intervals=source_cds,
        source_transcript_id=getattr(source_tx, "transcript_id", None),
    )

    kept_blocks, all_seqids = restrict_blocks_to_projected_transcript(projected_blocks, projected_tx)
    collapsed_cds = collapse_blocks_to_intervals(kept_blocks)

    seqid = projected_tx.seqid if projected_tx is not None else _dominant_seqid_from_blocks(kept_blocks)
    strand = projected_tx.strand if projected_tx is not None else _dominant_strand_from_blocks(kept_blocks, getattr(source_tx, "strand", None))

    projected_exons: list[tuple[int, int]] = []
    if projected_tx is not None:
        projected_exons = [normalize_interval(s, e) for s, e in projected_tx.exons]

    source_cds_seq = None
    if source_genome_by_seqid is not None and source_cds and getattr(source_tx, "seqid", None) in source_genome_by_seqid:
        source_cds_seq = extract_spliced_sequence(
            genome_by_seqid=source_genome_by_seqid,
            seqid=source_tx.seqid,
            intervals=source_cds,
            strand=source_tx.strand,
        )

    target_cds_seq = None
    if (
        target_genome_by_seqid is not None
        and seqid is not None
        and collapsed_cds
        and seqid in target_genome_by_seqid
    ):
        target_cds_seq = extract_spliced_sequence(
            genome_by_seqid=target_genome_by_seqid,
            seqid=seqid,
            intervals=collapsed_cds,
            strand=strand,
        )

    return ProjectedCDSModel(
        source_species=source_tx.species,
        source_transcript_id=getattr(source_tx, "transcript_id", getattr(source_tx, "source_transcript", "unknown")),
        target_species=target_species,
        target_gene_id=target_gene_id,
        target_transcript_id=target_transcript_id or getattr(projected_tx, "source_transcript", None),
        seqid=seqid,
        strand=strand,
        projected_exon_intervals=projected_exons,
        projected_cds_intervals=collapsed_cds,
        all_projected_cds_blocks=kept_blocks,
        projected_seqids=all_seqids,
        source_cds_length_nt=source_cds_len,
        target_cds_sequence=target_cds_seq,
        source_cds_sequence=source_cds_seq,
    )


# ---------------------------------------------------------------------------
# Assessment
# ---------------------------------------------------------------------------

def compute_basic_coding_metrics(model: ProjectedCDSModel) -> CodingIntegrityReport:
    seq = normalize_seq(model.target_cds_sequence)
    src_len = model.source_cds_length_nt
    phase_series = cds_phase_series(model.projected_cds_intervals, model.strand) if model.projected_cds_intervals else []

    if not seq:
        return CodingIntegrityReport(
            cds_length_nt=0,
            cds_length_aa=0,
            cds_recovery=0.0 if src_len else None,
            has_start_codon=None,
            has_stop_codon=None,
            internal_stop_count=None,
            length_mod_3=None,
            is_in_frame=None,
            phase_series=phase_series,
            source_cds_length_nt=src_len,
            target_cds_interval_count=len(model.projected_cds_intervals),
            classification_hint="missing_cds",
        )

    aa = translate_cds(seq)
    cds_len = len(seq)
    cds_recovery = None if src_len is None else cds_len / max(1, src_len)
    has_start = seq[:3] in START_CODONS if len(seq) >= 3 else False
    has_stop = seq[-3:] in STOP_CODONS if len(seq) >= 3 else False
    internal_stops = count_internal_stops(aa)
    mod3 = cds_len % 3
    in_frame = mod3 == 0

    if internal_stops > 0:
        hint = "internal_stops"
    elif not in_frame:
        hint = "out_of_frame"
    elif not has_start and not has_stop:
        hint = "missing_start_and_stop"
    elif not has_start:
        hint = "missing_start"
    elif not has_stop:
        hint = "missing_stop"
    elif cds_recovery is not None and cds_recovery < 0.5:
        hint = "low_recovery"
    else:
        hint = "apparently_intact"

    return CodingIntegrityReport(
        cds_length_nt=cds_len,
        cds_length_aa=len(aa),
        cds_recovery=cds_recovery,
        has_start_codon=has_start,
        has_stop_codon=has_stop,
        internal_stop_count=internal_stops,
        length_mod_3=mod3,
        is_in_frame=in_frame,
        phase_series=phase_series,
        source_cds_length_nt=src_len,
        target_cds_interval_count=len(model.projected_cds_intervals),
        classification_hint=hint,
    )


def detect_basic_fragmentation(model: ProjectedCDSModel) -> AssemblyFragmentationReport:
    seqids = sorted(set(model.projected_seqids))
    n_seqids = len(seqids)
    is_fragmented = n_seqids > 1
    hint = "multi_seqid_fragmentation" if is_fragmented else "single_seqid"
    return AssemblyFragmentationReport(
        is_fragmented_across_seqids=is_fragmented,
        n_seqids=n_seqids,
        seqids=seqids,
        classification_hint=hint,
    )


def classify_projected_locus(
    coding: CodingIntegrityReport,
    fragmentation: AssemblyFragmentationReport,
) -> tuple[str, str]:
    if coding.classification_hint == "missing_cds":
        return "lost_cds", "no projected CDS recovered"

    if coding.internal_stop_count and coding.internal_stop_count > 0:
        return "lost_cds", "internal stop codons detected"

    if coding.is_in_frame is False:
        return "lost_cds", "CDS length not divisible by 3"

    if fragmentation.is_fragmented_across_seqids:
        if coding.cds_recovery is not None and coding.cds_recovery >= 0.7:
            return "fragmented_cds", "split across multiple target seqids but substantial CDS recovered"
        return "fragmented_cds", "split across multiple target seqids"

    if (
        coding.cds_recovery is not None
        and coding.cds_recovery >= 0.95
        and coding.has_start_codon
        and coding.has_stop_codon
        and coding.internal_stop_count == 0
        and coding.is_in_frame
    ):
        return "intact_cds", "CDS appears intact"

    if (
        coding.cds_recovery is not None
        and coding.cds_recovery >= 0.7
        and coding.internal_stop_count == 0
        and coding.is_in_frame
    ):
        return "partial_cds", "substantial CDS recovered without obvious coding-disrupting defects"

    return "uncertain_cds", f"uncertain projection: {coding.classification_hint}"


def assess_projected_cds(model: ProjectedCDSModel) -> ProjectedLocusAssessment:
    coding = compute_basic_coding_metrics(model)
    fragmentation = detect_basic_fragmentation(model)
    final_class, final_reason = classify_projected_locus(coding, fragmentation)
    return ProjectedLocusAssessment(
        coding=coding,
        fragmentation=fragmentation,
        final_class=final_class,
        final_reason=final_reason,
    )


# ---------------------------------------------------------------------------
# GFF3 export helper
# ---------------------------------------------------------------------------

def gff3_attrs(attrs: dict[str, Any]) -> str:
    parts: list[str] = []
    for k, v in attrs.items():
        if v is None:
            continue
        if isinstance(v, bool):
            v = "true" if v else "false"
        elif isinstance(v, float):
            v = f"{v:.3f}"
        elif isinstance(v, list):
            v = ",".join(str(x) for x in v)
        else:
            v = str(v)
        parts.append(f"{k}={v}")
    return ";".join(parts)


def projected_cds_to_gff3_lines(
    model: ProjectedCDSModel,
    assessment: ProjectedLocusAssessment,
    source: str = "URCAT",
) -> list[str]:
    if model.seqid is None or not model.projected_cds_intervals:
        return []

    gene_id = model.target_gene_id or (model.target_transcript_id or model.source_transcript_id)
    tx_id = model.target_transcript_id or f"{gene_id}.1"

    transcript_span_intervals = model.projected_exon_intervals or model.projected_cds_intervals
    gene_start = min(s for s, _ in transcript_span_intervals)
    gene_end = max(e for _, e in transcript_span_intervals)

    attrs_common = {
        "source": source,
        "urcat_class": assessment.final_class,
        "urcat_reason": assessment.final_reason,
        "urcat_cds_recovery": assessment.coding.cds_recovery,
        "urcat_has_start_codon": assessment.coding.has_start_codon,
        "urcat_has_stop_codon": assessment.coding.has_stop_codon,
        "urcat_internal_stop_count": assessment.coding.internal_stop_count,
        "urcat_cds_length_nt": assessment.coding.cds_length_nt,
        "urcat_cds_length_aa": assessment.coding.cds_length_aa,
        "urcat_is_in_frame": assessment.coding.is_in_frame,
        "urcat_fragmented": assessment.fragmentation.is_fragmented_across_seqids,
        "urcat_fragment_seqids": assessment.fragmentation.seqids,
        "urcat_source_transcript": model.source_transcript_id,
    }

    lines: list[str] = []
    lines.append(
        "\t".join(
            [
                model.seqid,
                source,
                "gene",
                str(gene_start),
                str(gene_end),
                ".",
                model.strand or ".",
                ".",
                gff3_attrs({"ID": gene_id, "Name": gene_id, **attrs_common}),
            ]
        )
    )
    lines.append(
        "\t".join(
            [
                model.seqid,
                source,
                "mRNA",
                str(gene_start),
                str(gene_end),
                ".",
                model.strand or ".",
                ".",
                gff3_attrs({"ID": tx_id, "Parent": gene_id, "Name": tx_id, **attrs_common}),
            ]
        )
    )

    exon_intervals = model.projected_exon_intervals or model.projected_cds_intervals
    for i, (start, end) in enumerate(sorted(exon_intervals), start=1):
        lines.append(
            "\t".join(
                [
                    model.seqid,
                    source,
                    "exon",
                    str(start),
                    str(end),
                    ".",
                    model.strand or ".",
                    ".",
                    gff3_attrs({"ID": f"{tx_id}.exon{i}", "Parent": tx_id}),
                ]
            )
        )

    phases = cds_phase_series(model.projected_cds_intervals, model.strand)
    for i, ((start, end), phase) in enumerate(zip(sorted(model.projected_cds_intervals), phases), start=1):
        lines.append(
            "\t".join(
                [
                    model.seqid,
                    source,
                    "CDS",
                    str(start),
                    str(end),
                    ".",
                    model.strand or ".",
                    str(phase),
                    gff3_attrs({"ID": f"{tx_id}.cds{i}", "Parent": tx_id}),
                ]
            )
        )

    return lines
