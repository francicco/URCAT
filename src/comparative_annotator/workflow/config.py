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


def _require(cfg: configparser.ConfigParser, section: str, key: str) -> str:
    if not cfg.has_section(section):
        raise ValueError(f"Missing required section [{section}]")

    value = cfg.get(section, key, fallback="").strip()
    if not value:
        raise ValueError(f"Missing required key '{key}' in section [{section}]")
    return value


def _optional(cfg: configparser.ConfigParser, section: str, key: str, default: str = "") -> str:
    if not cfg.has_section(section):
        return default
    return cfg.get(section, key, fallback=default).strip()


def get_species_list_from_hal(hal_path: str) -> list[str]:
    """
    Parse species/genome names from halStats output.

    This expects halStats to print a line like:
        Genomes: Hmel Eisa Etal Diul
    """
    result = subprocess.run(
        ["halStats", hal_path],
        capture_output=True,
        text=True,
        check=True,
    )

    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("Genomes:"):
            genomes = line.split(":", 1)[1].strip()
            species = [x.strip() for x in genomes.split() if x.strip()]
            if not species:
                break
            return species

    raise RuntimeError(f"Could not parse species list from halStats output for: {hal_path}")


def _load_annotation_paths(
    cfg: configparser.ConfigParser,
    config_dir: Path,
) -> dict[str, str]:
    """
    Load [annotation] entries as explicit species -> GFF3 path mappings.

    Example:
        [annotation]
        Hmel = data/Hmel.test.gff3
        Eisa = data/Eisa.test.gff3
    """
    annotation_paths: dict[str, str] = {}

    if not cfg.has_section("annotation"):
        return annotation_paths

    for species, raw_path in cfg.items("annotation"):
        species = species.strip()
        raw_path = raw_path.strip()
        if not species or not raw_path:
            continue

        path = Path(raw_path)
        if not path.is_absolute():
            path = (config_dir / path).resolve()

        annotation_paths[species] = str(path)

    return annotation_paths


def _load_evidence(
    cfg: configparser.ConfigParser,
    config_dir: Path,
    species_list: list[str],
) -> dict[str, dict[str, str]]:
    """
    Load [evidence] into nested structure:
        evidence[species][evidence_type] = path

    Example input:
        [evidence]
        Hmel_bam = data/Hmel.bam
        Eisa_bam = data/Eisa.bam

    Result:
        {
            "Hmel": {"bam": "/abs/path/data/Hmel.bam"},
            "Eisa": {"bam": "/abs/path/data/Eisa.bam"},
            ...
        }
    """
    evidence: dict[str, dict[str, str]] = {sp: {} for sp in species_list}

    if not cfg.has_section("evidence"):
        return evidence

    for key, raw_path in cfg.items("evidence"):
        key = key.strip()
        raw_path = raw_path.strip()
        if not key or not raw_path:
            continue

        matched_species = None
        for sp in species_list:
            prefix = f"{sp}_"
            if key.startswith(prefix):
                matched_species = sp
                evidence_type = key[len(prefix):]
                break

        if matched_species is None or not evidence_type:
            continue

        path = Path(raw_path)
        if not path.is_absolute():
            path = (config_dir / path).resolve()

        evidence[matched_species][evidence_type] = str(path)

    return evidence


def load_urcat_config(config_path: str) -> URCATConfig:
    config_path_p = Path(config_path).resolve()
    config_dir = config_path_p.parent

    if not config_path_p.exists():
        raise FileNotFoundError(f"URCAT config file not found: {config_path_p}")

    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(config_path_p)

    seed_species = _require(cfg, "input", "seedSpecies")
    hal_path_raw = _require(cfg, "input", "halPath")
    batch_size_raw = _optional(cfg, "input", "batchSize", default="200")

    hal_path_p = Path(hal_path_raw)
    if not hal_path_p.is_absolute():
        hal_path_p = (config_dir / hal_path_p).resolve()

    try:
        batch_size = int(batch_size_raw)
    except ValueError as e:
        raise ValueError(f"Invalid integer for [input] batchSize: {batch_size_raw}") from e

    species_list = get_species_list_from_hal(str(hal_path_p))
    if seed_species not in species_list:
        raise ValueError(
            f"seedSpecies '{seed_species}' is not present in HAL species list: {species_list}"
        )

    annotation_paths = _load_annotation_paths(cfg, config_dir)
    evidence = _load_evidence(cfg, config_dir, species_list)

    return URCATConfig(
        seed_species=seed_species,
        hal_path=str(hal_path_p),
        batch_size=batch_size,
        species_list=species_list,
        annotation_paths=annotation_paths,
        evidence=evidence,
    )
