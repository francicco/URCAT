from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG"}


@dataclass
class ProjectedCDSModel:
    source_species: str
    source_transcript_id: str
    target_species: str
    target_gene_id: str
    target_transcript_id: str
    seqid: str | None
    strand: str | None
    exon_intervals: list[tuple[int, int]] = field(default_factory=list)
    cds_intervals: list[tuple[int, int]] = field(default_factory=list)
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
    classification_hint: str


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


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def normalize_seq(seq: str | None) -> str | None:
    if seq is None:
        return None
    return "".join(seq.split()).upper()


def sort_intervals(intervals: list[tuple[int, int]], strand: str | None = None) -> list[tuple[int, int]]:
    if strand == "-":
        return sorted(intervals, key=lambda x: (x[0], x[1]), reverse=False)
    return sorted(intervals, key=lambda x: (x[0], x[1]))


def cds_intervals_in_transcript_order(
    intervals: list[tuple[int, int]],
    strand: str | None,
) -> list[tuple[int, int]]:
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    if strand == "-":
        return list(reversed(intervals))
    return intervals


def interval_length(start: int, end: int) -> int:
    return end - start + 1


def sum_interval_lengths(intervals: list[tuple[int, int]]) -> int:
    return sum(interval_length(s, e) for s, e in intervals)


def extract_subsequence(
    genome_by_seqid: dict[str, str],
    seqid: str,
    start: int,
    end: int,
) -> str:
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
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa.append(codon_table.get(codon, "X"))
    return "".join(aa)


def count_internal_stops(aa: str) -> int:
    if not aa:
        return 0
    if aa.endswith("*"):
        aa = aa[:-1]
    return aa.count("*")


def cds_phase_series(intervals: list[tuple[int, int]], strand: str | None) -> list[int]:
    ordered = cds_intervals_in_transcript_order(intervals, strand)
    phases_in_tx_order = []
    consumed = 0
    for start, end in ordered:
        phase = (3 - (consumed % 3)) % 3
        if consumed == 0:
            phase = 0
        phases_in_tx_order.append(phase)
        consumed += interval_length(start, end)

    if strand == "-":
        phases_in_genomic_order = list(reversed(phases_in_tx_order))
    else:
        phases_in_genomic_order = phases_in_tx_order
    return phases_in_genomic_order


def get_transcript_cds_intervals(tx: Any) -> list[tuple[int, int]]:
    """
    Tries several common attribute names.
    """
    for attr in ("cds_parts", "cds_intervals", "cds", "CDS", "cds_blocks"):
        value = getattr(tx, attr, None)
        if value:
            out = []
            for item in value:
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    out.append((int(item[0]), int(item[1])))
            if out:
                return sorted(out, key=lambda x: (x[0], x[1]))
    return []


def get_transcript_exon_intervals(tx: Any) -> list[tuple[int, int]]:
    for attr in ("exons", "exon_intervals", "exon_blocks"):
        value = getattr(tx, attr, None)
        if value:
            out = []
            for item in value:
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    out.append((int(item[0]), int(item[1])))
            if out:
                return sorted(out, key=lambda x: (x[0], x[1]))
    return []


def project_intervals_with_hal(
    hal: Any,
    source_species: str,
    target_species: str,
    source_seqid: str,
    source_strand: str,
    intervals: list[tuple[int, int]],
    source_transcript_id: str | None = None,
) -> list[dict[str, Any]]:
    """
    Returns one flattened record per projected block.
    Expected HALAdapter API:
      hal.project_interval(
          source_species=...,
          target_species=...,
          seqid=...,
          start=...,
          end=...,
          strand=...,
          source_transcript=...
      )

    Each returned interval is expected to provide at least:
      seqid, start, end, strand
    """
    out = []
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
                {
                    "source_start": source_start,
                    "source_end": source_end,
                    "target_seqid": h["seqid"],
                    "target_start": int(h["start"]),
                    "target_end": int(h["end"]),
                    "target_strand": h.get("strand", source_strand),
                }
            )
    return out


def collapse_projected_blocks(
    projected_blocks: list[dict[str, Any]],
) -> tuple[str | None, str | None, list[tuple[int, int]], list[str]]:
    """
    Conservative rule:
    - if projected blocks land on multiple target seqids, keep all intervals but flag fragmentation
    - choose dominant seqid by total projected bp
    - choose dominant strand by support count
    """
    if not projected_blocks:
        return None, None, [], []

    seqid_bp: dict[str, int] = {}
    strand_n: dict[str, int] = {}

    for b in projected_blocks:
        seqid = b["target_seqid"]
        strand = b["target_strand"]
        bp = interval_length(b["target_start"], b["target_end"])
        seqid_bp[seqid] = seqid_bp.get(seqid, 0) + bp
        strand_n[strand] = strand_n.get(strand, 0) + 1

    dominant_seqid = sorted(seqid_bp.items(), key=lambda x: (-x[1], x[0]))[0][0]
    dominant_strand = sorted(strand_n.items(), key=lambda x: (-x[1], x[0]))[0][0]

    kept = [
        (b["target_start"], b["target_end"])
        for b in projected_blocks
        if b["target_seqid"] == dominant_seqid
    ]
    kept = sorted(kept, key=lambda x: (x[0], x[1]))

    all_seqids = sorted(seqid_bp.keys())
    return dominant_seqid, dominant_strand, kept, all_seqids


def build_projected_cds_model(
    source_tx: Any,
    target_species: str,
    target_gene_id: str,
    target_transcript_id: str,
    hal: Any,
    target_genome_by_seqid: dict[str, str] | None = None,
    source_genome_by_seqid: dict[str, str] | None = None,
) -> ProjectedCDSModel:
    source_cds = get_transcript_cds_intervals(source_tx)
    source_exons = get_transcript_exon_intervals(source_tx)

    projected_cds_blocks = project_intervals_with_hal(
        hal=hal,
        source_species=source_tx.species,
        target_species=target_species,
        source_seqid=source_tx.seqid,
        source_strand=source_tx.strand,
        intervals=source_cds,
        source_transcript_id=source_tx.transcript_id,
    )

    target_seqid, target_strand, target_cds_intervals, _ = collapse_projected_blocks(projected_cds_blocks)

    source_cds_seq = None
    target_cds_seq = None

    if source_genome_by_seqid is not None and source_cds and source_tx.seqid in source_genome_by_seqid:
        source_cds_seq = extract_spliced_sequence(
            genome_by_seqid=source_genome_by_seqid,
            seqid=source_tx.seqid,
            intervals=source_cds,
            strand=source_tx.strand,
        )

    if target_genome_by_seqid is not None and target_seqid is not None and target_cds_intervals:
        if target_seqid in target_genome_by_seqid:
            target_cds_seq = extract_spliced_sequence(
                genome_by_seqid=target_genome_by_seqid,
                seqid=target_seqid,
                intervals=target_cds_intervals,
                strand=target_strand,
            )

    return ProjectedCDSModel(
        source_species=source_tx.species,
        source_transcript_id=source_tx.transcript_id,
        target_species=target_species,
        target_gene_id=target_gene_id,
        target_transcript_id=target_transcript_id,
        seqid=target_seqid,
        strand=target_strand,
        exon_intervals=source_exons,
        cds_intervals=target_cds_intervals,
        source_cds_length_nt=sum_interval_lengths(source_cds) if source_cds else None,
        target_cds_sequence=target_cds_seq,
        source_cds_sequence=source_cds_seq,
    )


def compute_basic_coding_metrics(model: ProjectedCDSModel) -> CodingIntegrityReport:
    seq = normalize_seq(model.target_cds_sequence)
    src_len = model.source_cds_length_nt

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
            classification_hint="missing_cds",
        )

    aa = translate_cds(seq)
    cds_len = len(seq)
    cds_recovery = None if not src_len else cds_len / max(1, src_len)
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
        classification_hint=hint,
    )


def detect_basic_fragmentation(model: ProjectedCDSModel, all_projected_seqids: list[str] | None = None) -> AssemblyFragmentationReport:
    seqids = sorted(set(all_projected_seqids or ([] if model.seqid is None else [model.seqid])))
    n_seqids = len(seqids)
    is_fragmented = n_seqids > 1

    if is_fragmented:
        hint = "multi_seqid_fragmentation"
    else:
        hint = "single_seqid"

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
    """
    First-pass labels, intentionally simple.

    I   = intact
    PI  = partially intact
    UL  = uncertain/low-confidence
    L   = likely lost/pseudogenic
    AF  = assembly-fragmented
    """
    if fragmentation.is_fragmented_across_seqids:
        if coding.cds_recovery is not None and coding.cds_recovery >= 0.7:
            return "AF", "split across multiple target seqids but substantial CDS recovered"
        return "AF", "split across multiple target seqids"

    if coding.classification_hint == "missing_cds":
        return "L", "no projected CDS recovered"

    if coding.internal_stop_count and coding.internal_stop_count > 0:
        return "L", "internal stop codons detected"

    if coding.is_in_frame is False:
        return "L", "CDS length not divisible by 3"

    if (
        coding.cds_recovery is not None
        and coding.cds_recovery >= 0.95
        and coding.has_start_codon
        and coding.has_stop_codon
        and coding.internal_stop_count == 0
        and coding.is_in_frame
    ):
        return "I", "CDS appears intact"

    if (
        coding.cds_recovery is not None
        and coding.cds_recovery >= 0.7
        and coding.internal_stop_count == 0
        and coding.is_in_frame
    ):
        return "PI", "substantial CDS recovered without obvious coding-disrupting defects"

    return "UL", f"uncertain projection: {coding.classification_hint}"


def assess_projected_cds(
    model: ProjectedCDSModel,
    all_projected_seqids: list[str] | None = None,
) -> ProjectedLocusAssessment:
    coding = compute_basic_coding_metrics(model)
    fragmentation = detect_basic_fragmentation(model, all_projected_seqids=all_projected_seqids)
    final_class, final_reason = classify_projected_locus(coding, fragmentation)
    return ProjectedLocusAssessment(
        coding=coding,
        fragmentation=fragmentation,
        final_class=final_class,
        final_reason=final_reason,
    )


def gff3_attrs(attrs: dict[str, Any]) -> str:
    parts = []
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
    if model.seqid is None or not model.cds_intervals:
        return []

    gene_id = model.target_gene_id
    tx_id = model.target_transcript_id

    gene_start = min(s for s, _ in model.cds_intervals)
    gene_end = max(e for _, e in model.cds_intervals)

    attrs_common = {
        "source": source,
        "urcat_status": "new_locus",
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

    lines = []
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
                gff3_attrs(
                    {
                        "ID": gene_id,
                        "Name": gene_id,
                        **attrs_common,
                    }
                ),
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
                gff3_attrs(
                    {
                        "ID": tx_id,
                        "Parent": gene_id,
                        "Name": tx_id,
                        **attrs_common,
                    }
                ),
            ]
        )
    )

    phases = cds_phase_series(model.cds_intervals, model.strand)
    cds_intervals_sorted = sorted(model.cds_intervals, key=lambda x: (x[0], x[1]))

    for i, ((start, end), phase) in enumerate(zip(cds_intervals_sorted, phases), start=1):
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
                    gff3_attrs(
                        {
                            "ID": f"{tx_id}.cds{i}",
                            "Parent": tx_id,
                        }
                    ),
                ]
            )
        )

    return lines
