from __future__ import annotations

from pathlib import Path


def _read_noncomment_lines(path: str | Path) -> list[str]:
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            out.append(line)
    return out


def _feature_sort_key(line: str):
    parts = line.split("\t")
    seqid = parts[0]
    feature_type = parts[2]
    start = int(parts[3])
    end = int(parts[4])

    rank = {
        "gene": 0,
        "mRNA": 1,
        "transcript": 1,
        "exon": 2,
        "CDS": 3,
    }.get(feature_type, 9)

    return (seqid, start, end, rank, line)


def write_final_species_gff3s(
    workdir: str,
    annotation_dir: str,
    annotation_suffix: str,
    species_csv: str,
) -> None:
    species_list = [x.strip() for x in species_csv.split(",") if x.strip()]
    workdir = Path(workdir)
    out_dir = workdir / "final_gff3"
    out_dir.mkdir(parents=True, exist_ok=True)

    for species in species_list:
        original_path = Path(annotation_dir) / f"{species}{annotation_suffix}"
        species_new_files = sorted(workdir.rglob(f"{species}.new_loci.gff3"))

        original_lines = _read_noncomment_lines(original_path) if original_path.exists() else []
        new_lines = []
        for p in species_new_files:
            new_lines.extend(_read_noncomment_lines(p))

        merged_lines = sorted(original_lines + new_lines, key=_feature_sort_key)

        final_path = out_dir / f"{species}.final.gff3"
        new_only_path = out_dir / f"{species}.new_loci.gff3"

        with open(final_path, "w") as fh:
            fh.write("##gff-version 3\n")
            for line in merged_lines:
                fh.write(line + "\n")

        with open(new_only_path, "w") as fh:
            fh.write("##gff-version 3\n")
            for line in sorted(new_lines, key=_feature_sort_key):
                fh.write(line + "\n")
