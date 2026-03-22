from __future__ import annotations

import subprocess
import uuid
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
    out_fasta: str,
) -> str:
    """
    Export one species genome from HAL to FASTA.
    """
    out = Path(out_fasta)
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


def run_gffread(gff_path: str, genome_fa: str, prefix: str) -> tuple[str, str, str]:
    """
    Extract mRNA, CDS, and protein FASTA with gffread.

    Some runs may fail to emit one or more outputs for edge cases. To keep the
    pipeline robust, missing outputs are replaced by empty files.
    """
    prefix_path = Path(prefix)
    _safe_mkdir(prefix_path.parent)

    mrna = Path(f"{prefix}.mrna.fa")
    cds = Path(f"{prefix}.cds.fa")
    aa = Path(f"{prefix}.aa.fa")

    tmp_tag = uuid.uuid4().hex
    mrna_tmp = Path(f"{prefix}.mrna.fa.{tmp_tag}.tmp")
    cds_tmp = Path(f"{prefix}.cds.fa.{tmp_tag}.tmp")
    aa_tmp = Path(f"{prefix}.aa.fa.{tmp_tag}.tmp")

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

    for tmp in (mrna_tmp, cds_tmp, aa_tmp):
        if not tmp.exists():
            tmp.write_text("")

    mrna_tmp.replace(mrna)
    cds_tmp.replace(cds)
    aa_tmp.replace(aa)

    return str(mrna), str(cds), str(aa)


def sanitize_protein_fasta(in_fa: str, out_fa: str | None = None) -> str:
    """
    Replace DIAMOND-invalid protein symbols with X.

    In practice this mainly fixes '.' characters, but it also sanitizes any
    other unexpected residue symbols.
    """
    in_path = Path(in_fa)
    out_path = Path(out_fa) if out_fa is not None else in_path
    tmp_path = out_path.with_suffix(out_path.suffix + ".sanitize.tmp")

    valid = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

    with open(in_path) as fin, open(tmp_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                seq = line.strip().upper().replace(".", "X")
                seq = "".join(ch if ch in valid else "X" for ch in seq)
                fout.write(seq + "\n")

    tmp_path.replace(out_path)
    return str(out_path)


def prepare_species_sequences(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_list: list[str],
    hal_path: str,
) -> dict[str, dict[str, str | bool | None]]:
    """
    For each species:
    - export genome FASTA from HAL
    - if annotation exists, derive mRNA/CDS/AA FASTA with gffread
    - sanitize AA FASTA for DIAMOND
    """
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

        sanitize_protein_fasta(aa_fa)

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
    cache_dir: str,
    source_species: str,
    target_species: str,
    sequences_by_species: dict[str, dict[str, str | bool | None]],
) -> tuple[str, str]:
    """
    Return (query_fa, target_fa) for one source->target DIAMOND run.
    """
    source = sequences_by_species.get(source_species)
    target = sequences_by_species.get(target_species)

    if source is None:
        raise ValueError(f"Missing sequences for source species: {source_species}")
    if target is None:
        raise ValueError(f"Missing sequences for target species: {target_species}")

    query_fa = source.get("aa_fa")
    target_fa = target.get("aa_fa")

    if not query_fa:
        raise ValueError(f"Source species has no protein FASTA: {source_species}")
    if not target_fa:
        raise ValueError(f"Target species has no protein FASTA: {target_species}")

    return str(query_fa), str(target_fa)


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
    """
    out_tsv_p = Path(out_tsv)
    db_prefix = Path(f"{tmp_prefix}.dmnd")

    _safe_mkdir(out_tsv_p.parent)
    _safe_mkdir(db_prefix.parent)

    sanitize_protein_fasta(target_fa)
    sanitize_protein_fasta(query_fa)

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
    Read DIAMOND outfmt 6 into:
      {(qseqid, sseqid): {"pid": ..., "aln_len": ..., "bitscore": ..., "evalue": ...}}
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
