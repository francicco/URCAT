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

def sanitize_protein_sequence(seq: str) -> tuple[str, dict]:
    seq = seq.strip().upper()

    flags = {
        "had_terminal_stop": False,
        "had_internal_stop": False,
        "had_invalid_chars": False,
        "is_empty_after_cleaning": False,
    }

    if seq.endswith(".") or seq.endswith("*"):
        flags["had_terminal_stop"] = True
        seq = seq[:-1]

    if "." in seq or "*" in seq:
        flags["had_internal_stop"] = True
        seq = seq.replace(".", "").replace("*", "")

    allowed = set("ACDEFGHIKLMNPQRSTVWYX")
    cleaned = []
    for aa in seq:
        if aa in allowed:
            cleaned.append(aa)
        else:
            flags["had_invalid_chars"] = True
            cleaned.append("X")

    seq = "".join(cleaned)

    if not seq:
        flags["is_empty_after_cleaning"] = True

    return seq, flags

def read_protein_fasta_with_qc(path: str) -> tuple[dict[str, str], dict[str, dict]]:
    seqs = {}
    qc = {}

    name = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raw = "".join(chunks)
                    clean, flags = sanitize_protein_sequence(raw)
                    seqs[name] = clean
                    qc[name] = flags
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        raw = "".join(chunks)
        clean, flags = sanitize_protein_sequence(raw)
        seqs[name] = clean
        qc[name] = flags

    return seqs, qc

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
            aa_seqs, aa_qc = read_protein_fasta_with_qc(paths["aa"]),
        }

    return seqs

def write_fasta(path: str, seqs: dict[str, str]):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            if not seq:
                continue
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def filter_aa_for_diamond(aa_seqs: dict[str, str], aa_qc: dict[str, dict]) -> dict[str, str]:
    out = {}
    for tx_id, seq in aa_seqs.items():
        flags = aa_qc.get(tx_id, {})
        if flags.get("is_empty_after_cleaning"):
            continue
        if flags.get("had_internal_stop"):
            continue
        out[tx_id] = seq
    return out
