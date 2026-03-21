from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from comparative_annotator.workflow.progressive import read_json


@dataclass
class ConsensusLocus:
    species: str
    seqid: str
    strand: str
    support_count: int
    total_chain_score: float
    mean_exon_recovery: float
    source_transcripts: list[str]
    exons: list[tuple[int, int]]


def _consensus_from_dict(d: dict[str, Any]) -> ConsensusLocus:
    return ConsensusLocus(
        species=d["species"],
        seqid=d["seqid"],
        strand=d["strand"],
        support_count=int(d.get("support_count", 0)),
        total_chain_score=float(d.get("total_chain_score", 0.0)),
        mean_exon_recovery=float(d.get("mean_exon_recovery", 0.0)),
        source_transcripts=list(d.get("source_transcripts", [])),
        exons=[(int(a), int(b)) for a, b in d.get("exons", [])],
    )


def _sort_key(c: ConsensusLocus):
    if c.exons:
        start = min(a for a, _ in c.exons)
        end = max(b for _, b in c.exons)
    else:
        start = 10**18
        end = 10**18
    return (c.seqid, start, end, c.strand, tuple(c.exons), tuple(sorted(c.source_transcripts)))


def _dedup_key(c: ConsensusLocus):
    return (
        c.species,
        c.seqid,
        c.strand,
        tuple(c.exons),
    )


def deduplicate_consensus_loci(consensuses: list[ConsensusLocus]) -> list[ConsensusLocus]:
    best_by_key: dict[tuple, ConsensusLocus] = {}

    for c in consensuses:
        key = _dedup_key(c)
        prev = best_by_key.get(key)

        if prev is None:
            best_by_key[key] = c
            continue

        prev_score = (prev.support_count, prev.total_chain_score, prev.mean_exon_recovery)
        cur_score = (c.support_count, c.total_chain_score, c.mean_exon_recovery)

        if cur_score > prev_score:
            best_by_key[key] = c

    out = list(best_by_key.values())
    out.sort(key=_sort_key)
    return out


def write_new_loci_gff3(output_path: str | Path, consensuses: list[ConsensusLocus]) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    consensuses = [c for c in consensuses if c.exons]
    consensuses = deduplicate_consensus_loci(consensuses)

    with open(output_path, "w") as out:
        out.write("##gff-version 3\n")

        for i, c in enumerate(consensuses, start=1):
            gene_id = f"{c.species}.URCAT.newG{i:06d}"
            tx_id = f"{gene_id}.1"

            start = min(a for a, _ in c.exons)
            end = max(b for _, b in c.exons)

            gene_attrs = {
                "ID": gene_id,
                "Name": gene_id,
                "source": "URCAT",
                "urcat_status": "new_locus",
                "urcat_support_count": str(c.support_count),
                "urcat_total_chain_score": f"{c.total_chain_score:.3f}",
                "urcat_mean_exon_recovery": f"{c.mean_exon_recovery:.3f}",
                "urcat_source_transcripts": ",".join(c.source_transcripts),
            }
            mrna_attrs = {
                "ID": tx_id,
                "Parent": gene_id,
                "Name": tx_id,
                "source": "URCAT",
                "urcat_status": "new_locus",
                "urcat_support_count": str(c.support_count),
                "urcat_total_chain_score": f"{c.total_chain_score:.3f}",
                "urcat_mean_exon_recovery": f"{c.mean_exon_recovery:.3f}",
                "urcat_source_transcripts": ",".join(c.source_transcripts),
            }

            out.write(
                f"{c.seqid}\tURCAT\tgene\t{start}\t{end}\t.\t{c.strand}\t.\t"
                + ";".join(f"{k}={v}" for k, v in gene_attrs.items())
                + "\n"
            )
            out.write(
                f"{c.seqid}\tURCAT\tmRNA\t{start}\t{end}\t.\t{c.strand}\t.\t"
                + ";".join(f"{k}={v}" for k, v in mrna_attrs.items())
                + "\n"
            )

            for exon_i, (exon_start, exon_end) in enumerate(c.exons, start=1):
                exon_id = f"{tx_id}.exon{exon_i}"
                exon_attrs = {
                    "ID": exon_id,
                    "Parent": tx_id,
                }
                out.write(
                    f"{c.seqid}\tURCAT\texon\t{exon_start}\t{exon_end}\t.\t{c.strand}\t.\t"
                    + ";".join(f"{k}={v}" for k, v in exon_attrs.items())
                    + "\n"
                )

    return str(output_path)


def write_round_new_loci_gff3(round_ref_dir: str | Path, decision_path: str | Path) -> dict[str, str]:
    round_ref_dir = Path(round_ref_dir)
    decision = read_json(decision_path)

    outputs: dict[str, str] = {}

    for species, items in sorted((decision.get("new_consensus_by_species") or {}).items()):
        consensuses = [_consensus_from_dict(x) for x in items]
        outpath = round_ref_dir / f"target_{species}" / f"{species}.new_loci.gff3"
        outputs[species] = write_new_loci_gff3(outpath, consensuses)

    return outputs


def collect_all_new_consensus_by_species(output_dir: str | Path) -> dict[str, list[ConsensusLocus]]:
    output_dir = Path(output_dir)
    by_species: dict[str, list[ConsensusLocus]] = {}

    for decision_path in sorted(output_dir.glob("rounds/round_*/ref_*/post_round_decision.json")):
        decision = read_json(decision_path)
        for species, items in (decision.get("new_consensus_by_species") or {}).items():
            by_species.setdefault(species, []).extend(_consensus_from_dict(x) for x in items)

    for species in list(by_species):
        by_species[species] = deduplicate_consensus_loci(by_species[species])

    return by_species


def copy_original_plus_new_to_final(
    output_dir: str | Path,
    annotation_dir: str | Path,
    annotation_suffix: str,
    species_list: list[str],
) -> dict[str, dict[str, str]]:
    output_dir = Path(output_dir)
    annotation_dir = Path(annotation_dir)
    final_dir = output_dir / "final_gff3"
    final_dir.mkdir(parents=True, exist_ok=True)

    consensus_by_species = collect_all_new_consensus_by_species(output_dir)
    outputs: dict[str, dict[str, str]] = {}

    for species in species_list:
        original_gff = annotation_dir / f"{species}{annotation_suffix}"
        new_only_gff = final_dir / f"{species}.new_loci.gff3"
        final_gff = final_dir / f"{species}.final.gff3"

        consensuses = consensus_by_species.get(species, [])
        write_new_loci_gff3(new_only_gff, consensuses)

        with open(final_gff, "w") as out:
            wrote_header = False

            if original_gff.exists():
                with open(original_gff) as fh:
                    for line in fh:
                        if line.startswith("##gff-version 3"):
                            if not wrote_header:
                                out.write(line)
                                wrote_header = True
                            continue
                        out.write(line)

            if not wrote_header:
                out.write("##gff-version 3\n")

            if new_only_gff.exists():
                with open(new_only_gff) as fh:
                    for line in fh:
                        if line.startswith("##gff-version 3"):
                            continue
                        out.write(line)

        outputs[species] = {
            "new_loci_gff3": str(new_only_gff),
            "final_gff3": str(final_gff),
        }

    return outputs
