from __future__ import annotations

import configparser
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class URCATConfig:
    seed_species: str
    hal_path: str
    batch_size: int
    species_list: list[str]
    annotation_paths: dict[str, str] = field(default_factory=dict)
    evidence: dict[str, dict[str, str]] = field(default_factory=dict)


def _require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def get_hal_species_list(hal_path: str) -> list[str]:
    """
    Return all species/genomes in the HAL file using halStats --genomes.
    """
    result = subprocess.run(
        ["halStats", "--genomes", str(hal_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    species = [x.strip() for x in result.stdout.strip().split() if x.strip()]
    if not species:
        raise RuntimeError(f"No species found in HAL: {hal_path}")
    return species


def _normalize_path(value: str, config_dir: Path) -> str:
    p = Path(value)
    if not p.is_absolute():
        p = (config_dir / p).resolve()
    else:
        p = p.resolve()
    return str(p)


def _parse_annotation_section(
    cp: configparser.ConfigParser,
    config_dir: Path,
    species_list: list[str],
) -> dict[str, str]:
    """
    Parse:
      [annotation]
      Hmel = data/Hmel.test.gff3
      Eisa = data/Eisa.test.gff3
    Missing species are allowed.
    """
    paths: dict[str, str] = {}

    if not cp.has_section("annotation"):
        return paths

    for species, raw_path in cp.items("annotation"):
        species_name = species.strip()
        if species_name not in species_list:
            continue
        paths[species_name] = _normalize_path(raw_path.strip(), config_dir)

    return paths


def _parse_evidence_section(
    cp: configparser.ConfigParser,
    config_dir: Path,
    species_list: list[str],
) -> dict[str, dict[str, str]]:
    """
    Parse:
      [evidence]
      Hmel_bam = data/Hmel.bam
      Eisa_bam = data/Eisa.bam

    Returns:
      {
        "Hmel": {"bam": "/abs/path/Hmel.bam"},
        "Eisa": {"bam": "/abs/path/Eisa.bam"},
      }
    """
    evidence: dict[str, dict[str, str]] = {sp: {} for sp in species_list}

    if not cp.has_section("evidence"):
        return evidence

    for key, raw_path in cp.items("evidence"):
        key = key.strip()
        value = raw_path.strip()

        if "_" not in key:
            continue

        species, evtype = key.split("_", 1)
        if species not in evidence:
            continue

        evidence[species][evtype] = _normalize_path(value, config_dir)

    return evidence


def load_urcat_config(config_path: str) -> URCATConfig:
    """
    Expected format:

    [input]
    seedSpecies = Hmel
    halPath = data/3SpChr21.hal
    batchSize = 5

    [annotation]
    Hmel = data/Hmel.test.gff3
    Eisa = data/Eisa.test.gff3
    Diul = data/Diul.test.gff3

    [evidence]
    Hmel_bam = data/Hmel.bam
    Eisa_bam = data/Eisa.bam
    Diul_bam = data/Diul.bam
    """
    config_file = Path(config_path).resolve()
    _require_file(config_file, "Config file")

    cp = configparser.ConfigParser()
    cp.optionxform = str
    cp.read(config_file)

    if not cp.has_section("input"):
        raise ValueError("Config file must contain an [input] section")

    input_section = cp["input"]

    seed_species = input_section.get("seedSpecies", "").strip()
    hal_path_raw = input_section.get("halPath", "").strip()
    batch_size_raw = input_section.get("batchSize", "200").strip()

    if not seed_species:
        raise ValueError("Missing required key [input] seedSpecies")
    if not hal_path_raw:
        raise ValueError("Missing required key [input] halPath")

    hal_path = _normalize_path(hal_path_raw, config_file.parent)
    _require_file(Path(hal_path), "HAL file")

    try:
        batch_size = int(batch_size_raw)
    except ValueError as e:
        raise ValueError(f"Invalid batchSize: {batch_size_raw}") from e

    species_list = get_hal_species_list(hal_path)

    if seed_species not in species_list:
        raise ValueError(
            f"seedSpecies '{seed_species}' is not present in HAL species list: {species_list}"
        )

    annotation_paths = _parse_annotation_section(cp, config_file.parent, species_list)
    evidence = _parse_evidence_section(cp, config_file.parent, species_list)

    return URCATConfig(
        seed_species=seed_species,
        hal_path=hal_path,
        batch_size=batch_size,
        species_list=species_list,
        annotation_paths=annotation_paths,
        evidence=evidence,
    )
