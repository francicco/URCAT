from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from comparative_annotator.models.transcript import CandidateTranscript


@dataclass
class GFF3Feature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: Dict[str, str]


def parse_attributes(attr: str) -> Dict[str, str]:
    """
    Parse a GFF3 attribute column into a dictionary.

    Example:
        ID=tx1;Parent=gene1;Name=ABC
    """
    out: Dict[str, str] = {}

    attr = attr.strip()
    if not attr:
        return out

    for item in attr.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        out[key] = value

    return out


def parse_gff3_line(line: str) -> Optional[GFF3Feature]:
    """
    Parse one non-comment GFF3 line.
    Returns None for malformed lines.
    """
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 9:
        return None

    seqid, source, feature_type, start, end, score, strand, phase, attr = parts

    try:
        start_i = int(start)
        end_i = int(end)
    except ValueError:
        return None

    return GFF3Feature(
        seqid=seqid,
        source=source,
        feature_type=feature_type,
        start=start_i,
        end=end_i,
        score=score,
        strand=strand,
        phase=phase,
        attributes=parse_attributes(attr),
    )


def _as_parent_list(attributes: Dict[str, str]) -> List[str]:
    parent = attributes.get("Parent")
    if not parent:
        return []
    return [x.strip() for x in parent.split(",") if x.strip()]


def load_gff3(
    path: str | Path,
    species: str,
    source_name: Optional[str] = None,
) -> Dict[str, CandidateTranscript]:
    """
    Load a GFF3 file and return transcript models keyed by transcript_id.

    Supported feature types:
      - gene
      - mRNA / transcript
      - exon
      - CDS

    Important:
      - exons/CDS may appear before transcript lines
      - exon/CDS features may have multiple parents
    """
    path = Path(path)
    source_name = source_name or path.stem

    # Transcript metadata, keyed by transcript ID
    transcript_meta: Dict[str, Dict[str, object]] = {}

    # Child features temporarily grouped by transcript ID
    exon_map: Dict[str, List[tuple[int, int]]] = {}
    cds_map: Dict[str, List[tuple[int, int]]] = {}

    # Optional gene parent lookup
    transcript_to_gene: Dict[str, str] = {}

    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            feature = parse_gff3_line(raw)
            if feature is None:
                continue

            ftype = feature.feature_type.lower()
            attrs = feature.attributes

            if ftype in {"mrna", "transcript"}:
                tx_id = attrs.get("ID")
                if not tx_id:
                    continue

                transcript_meta[tx_id] = {
                    "seqid": feature.seqid,
                    "start": feature.start,
                    "end": feature.end,
                    "strand": feature.strand,
                    "source": source_name,
                }

                parent_gene = attrs.get("Parent")
                if parent_gene:
                    transcript_to_gene[tx_id] = parent_gene

            elif ftype == "exon":
                parents = _as_parent_list(attrs)
                if not parents:
                    continue

                for tx_id in parents:
                    exon_map.setdefault(tx_id, []).append((feature.start, feature.end))

                    # If transcript line has not appeared yet, create placeholder metadata
                    if tx_id not in transcript_meta:
                        transcript_meta[tx_id] = {
                            "seqid": feature.seqid,
                            "start": feature.start,
                            "end": feature.end,
                            "strand": feature.strand,
                            "source": source_name,
                        }
                    else:
                        transcript_meta[tx_id]["start"] = min(
                            int(transcript_meta[tx_id]["start"]), feature.start
                        )
                        transcript_meta[tx_id]["end"] = max(
                            int(transcript_meta[tx_id]["end"]), feature.end
                        )

            elif ftype == "cds":
                parents = _as_parent_list(attrs)
                if not parents:
                    continue

                for tx_id in parents:
                    cds_map.setdefault(tx_id, []).append((feature.start, feature.end))

                    if tx_id not in transcript_meta:
                        transcript_meta[tx_id] = {
                            "seqid": feature.seqid,
                            "start": feature.start,
                            "end": feature.end,
                            "strand": feature.strand,
                            "source": source_name,
                        }
                    else:
                        transcript_meta[tx_id]["start"] = min(
                            int(transcript_meta[tx_id]["start"]), feature.start
                        )
                        transcript_meta[tx_id]["end"] = max(
                            int(transcript_meta[tx_id]["end"]), feature.end
                        )

            # We ignore gene lines for now, except through transcript Parent links.

    transcripts: Dict[str, CandidateTranscript] = {}

    for tx_id, meta in transcript_meta.items():
        tx = CandidateTranscript(
            transcript_id=tx_id,
            species=species,
            source=source_name,
            seqid=str(meta["seqid"]),
            start=int(meta["start"]),
            end=int(meta["end"]),
            strand=str(meta["strand"]),
            exons=exon_map.get(tx_id, []),
            cds=cds_map.get(tx_id, []),
        )

        tx.attributes["gene_id"] = transcript_to_gene.get(tx_id)
        tx.finalize()
        transcripts[tx_id] = tx

    return transcripts
