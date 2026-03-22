from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class URCATConfig:
    raw: dict[str, Any]

    @property
    def output_dir(self) -> str:
        return self.raw["project"]["output_dir"]

    @property
    def hal_path(self) -> str:
        return self.raw["project"]["hal_path"]

    @property
    def annotation_dir(self) -> str:
        return self.raw["project"]["annotation_dir"]

    @property
    def annotation_suffix(self) -> str:
        return self.raw["project"].get("annotation_suffix", ".gff3")

    @property
    def seed_species(self) -> str:
        return self.raw["species"]["seed_species"]

    @property
    def species_list(self) -> list[str]:
        return list(self.raw["species"]["all"])

    @property
    def annotation_map(self) -> dict[str, str | None]:
        return dict(self.raw["species"].get("annotations", {}))

    @property
    def batch_size(self) -> int:
        return int(self.raw["workflow"].get("batch_size", 200))

    @property
    def allow_missing_annotations(self) -> bool:
        return bool(self.raw["workflow"].get("allow_missing_annotations", True))


def load_config(path: str) -> URCATConfig:
    with open(path, "r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)

    if raw is None:
        raise ValueError("Empty config file")

    required_top = ["project", "species", "workflow"]
    for key in required_top:
        if key not in raw:
            raise ValueError(f"Missing top-level config section: {key}")

    cfg = URCATConfig(raw)

    if not Path(cfg.hal_path).exists():
        raise FileNotFoundError(f"HAL file not found: {cfg.hal_path}")

    if cfg.seed_species not in cfg.species_list:
        raise ValueError(f"seed_species {cfg.seed_species!r} is not in species.all")

    return cfg
