from __future__ import annotations

import hashlib
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class SpeciesSequencePaths:
    species: str
    genome_fa: str
    mrna_fa: str | None
    cds_fa: str | None
    aa_fa: str | None
    has_annotation: bool


def _safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _hash_file_identity(path: Path) -> str:
    return hashlib.md5(str(path.resolve()).encode("utf-8")).hexdigest()


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def ensure_fasta_index(fasta_path: str) -> None:
    fasta = Path(fasta_path)
    fai = Path(str(fasta) + ".fai")
    if fai.exists():
        return
    _run(["samtools", "faidx", str(fasta)])


def export_species_genome_from_hal(
    hal_path: str,
    species: str,
    out_fa: str,
) -> str:
    out = Path(out_fa)
    _safe_mkdir(out.parent)

    if out.exists() and out.stat().st_size > 0:
        return str(out)

    _run(
        [
            "hal2fasta",
            hal_path,
            species,
            "--outFaPath",
            str(out),
        ]
    )

    return str(out)

def run_gffread(
    gff_path: str,
    genome_fa: str,
    prefix: str,
) -> tuple[str, str, str]:
    prefix_path = Path(prefix)
    _safe_mkdir(prefix_path.parent)

    token = _hash_file_identity(Path(gff_path))
    mrna_tmp = prefix_path.parent / f"{prefix_path.name}.mrna.fa.{token}.tmp"
    cds_tmp = prefix_path.parent / f"{prefix_path.name}.cds.fa.{token}.tmp"
    aa_tmp = prefix_path.parent / f"{prefix_path.name}.aa.fa.{token}.tmp"

    mrna = prefix_path.parent / f"{prefix_path.name}.mrna.fa"
    cds = prefix_path.parent / f"{prefix_path.name}.cds.fa"
    aa = prefix_path.parent / f"{prefix_path.name}.aa.fa"

    if mrna.exists() and cds.exists() and aa.exists():
        return str(mrna), str(cds), str(aa)

    cmd = [
        "gffread",
        gff_path,
        "-g",
        genome_fa,
        "-w",
        str(mrna_tmp),
        "-x",
        str(cds_tmp),
        "-y",
        str(aa_tmp),
    ]
    _run(cmd)

    mrna_tmp.replace(mrna)
    cds_tmp.replace(cds)
    aa_tmp.replace(aa)

    return str(mrna), str(cds), str(aa)


def prepare_species_sequences(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_list: list[str],
    hal_path: str,
) -> dict[str, dict[str, str | bool | None]]:
    workdir_p = Path(workdir)
    annotation_dir_p = Path(annotation_dir)
    cache_dir = workdir_p / "sequence_cache"
    _safe_mkdir(cache_dir)

    out: dict[str, dict[str, str | bool | None]] = {}

    for species in species_list:
        genome_fa = cache_dir / f"{species}.fa"
        export_species_genome_from_hal(
            hal_path=hal_path,
            species=species,
            out_fasta=str(genome_fa),
        )

        gff = annotation_dir_p / f"{species}{annotation_suffix}"

        if not gff.exists():
            out[species] = {
                "genome_fa": str(genome_fa),
                "mrna_fa": None,
                "cds_fa": None,
                "aa_fa": None,
                "has_annotation": False,
            }
            continue

        prefix = cache_dir / species
        mrna_fa, cds_fa, aa_fa = run_gffread(
            gff_path=str(gff),
            genome_fa=str(genome_fa),
            prefix=str(prefix),
        )

        out[species] = {
            "genome_fa": str(genome_fa),
            "mrna_fa": mrna_fa,
            "cds_fa": cds_fa,
            "aa_fa": aa_fa,
            "has_annotation": True,
        }

    return out


def load_all_species_sequences(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_list: list[str],
    hal_path: str,
) -> dict[str, dict[str, str | bool | None]]:
    return prepare_species_sequences(
        workdir=workdir,
        annotation_dir=annotation_dir,
        annotation_suffix=annotation_suffix,
        species_list=species_list,
        hal_path=hal_path,
    )


def prepare_diamond_inputs(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_list: list[str],
    hal_path: str,
) -> dict[str, str]:
    """
    Return only protein FASTA paths for annotated species.
    Unannotated species are skipped.
    """
    seqs = load_all_species_sequences(
        workdir=workdir,
        annotation_dir=annotation_dir,
        annotation_suffix=annotation_suffix,
        species_list=species_list,
        hal_path=hal_path,
    )

    return {
        species: paths["aa_fa"]
        for species, paths in seqs.items()
        if paths.get("aa_fa") is not None
    }


def run_diamond(
    query_fa: str,
    target_fa: str,
    out_tsv: str,
    tmp_prefix: str,
    sensitivity: str = "--ultra-sensitive",
    max_target_seqs: int = 25,
    evalue: float = 1e-5,
    threads: int = 1,
) -> str:
    """
    Run DIAMOND blastp of query proteins against target proteins.
    Returns the output TSV path.
    """
    out_tsv_p = Path(out_tsv)
    db_prefix = Path(f"{tmp_prefix}.dmnd")

    _safe_mkdir(out_tsv_p.parent)
    _safe_mkdir(db_prefix.parent)

    if not db_prefix.exists():
        _run(
            [
                "diamond",
                "makedb",
                "--in",
                str(target_fa),
                "--db",
                str(db_prefix),
            ]
        )

    cmd = [
        "diamond",
        "blastp",
        "--query",
        str(query_fa),
        "--db",
        str(db_prefix),
        "--out",
        str(out_tsv_p),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "bitscore",
        "evalue",
        "--threads",
        str(threads),
        "--max-target-seqs",
        str(max_target_seqs),
        "--evalue",
        str(evalue),
    ]

    if sensitivity:
        cmd.insert(2, sensitivity)

    _run(cmd)
    return str(out_tsv_p)


def load_diamond_results(path: str) -> dict[tuple[str, str], dict]:
    """
    Read DIAMOND outfmt 6 table into:
      {(qseqid, sseqid): {"pid": ..., "aln_len": ..., "bitscore": ..., "evalue": ...}}
    keeping the first/best row encountered per pair.
    """
    results: dict[tuple[str, str], dict] = {}
    p = Path(path)

    if not p.exists() or p.stat().st_size == 0:
        return results

    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 6:
                continue

            qseqid, sseqid, pident, length, bitscore, evalue = fields[:6]
            key = (qseqid, sseqid)

            if key not in results:
                results[key] = {
                    "pid": float(pident),
                    "aln_len": int(length),
                    "bitscore": float(bitscore),
                    "evalue": float(evalue),
                }

    return results
