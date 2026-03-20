from __future__ import annotations

import subprocess
from pathlib import Path


def read_fasta(path: str) -> dict[str, str]:
    seqs = {}
    name = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        seqs[name] = "".join(chunks)

    return seqs


def hal_to_fasta(hal_path: str, species: str, out_fa: str):
    if Path(out_fa).exists():
        return

    cmd = [
        "hal2fasta",
        hal_path,
        species,
    ]

    with open(out_fa, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def run_gffread(gff: str, genome_fa: str, prefix: str):
    mrna = f"{prefix}.mrna.fa"
    cds = f"{prefix}.cds.fa"
    aa = f"{prefix}.aa.fa"

    if Path(aa).exists() and Path(cds).exists() and Path(mrna).exists():
        return mrna, cds, aa

    cmd = [
        "gffread",
        gff,
        "-g", genome_fa,
        "-w", mrna,
        "-x", cds,
        "-y", aa,
    ]
    subprocess.run(cmd, check=True)

    return mrna, cds, aa


def prepare_species_sequences(
    hal_path: str,
    annotation_dir: str,
    annotation_suffix: str,
    cache_dir: str,
    species: str,
):
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    genome_fa = cache_dir / f"{species}.fa"
    hal_to_fasta(hal_path, species, str(genome_fa))

    gff = Path(annotation_dir) / f"{species}{annotation_suffix}"

    prefix = str(cache_dir / species)
    mrna, cds, aa = run_gffread(str(gff), str(genome_fa), prefix)

    return {
        "mrna": mrna,
        "cds": cds,
        "aa": aa,
    }


def load_all_species_sequences(
    hal_path: str,
    annotation_dir: str,
    annotation_suffix: str,
    cache_dir: str,
    species_list: list[str],
):
    seqs = {}

    for sp in species_list:
        paths = prepare_species_sequences(
            hal_path=hal_path,
            annotation_dir=annotation_dir,
            annotation_suffix=annotation_suffix,
            cache_dir=cache_dir,
            species=sp,
        )

        seqs[sp] = {
            "mrna": read_fasta(paths["mrna"]),
            "cds": read_fasta(paths["cds"]),
            "aa": read_fasta(paths["aa"]),
        }

    return seqs
