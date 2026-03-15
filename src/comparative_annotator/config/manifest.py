from dataclasses import dataclass
from pathlib import Path
import yaml


@dataclass
class SpeciesConfig:
    name: str
    genome: Path | None = None
    annotations: list[Path] | None = None
    proteins: list[Path] | None = None


@dataclass
class ProjectConfig:
    hal: Path
    species: dict[str, SpeciesConfig]
    references: list[str]


def load_manifest(path: str | Path) -> ProjectConfig:

    with open(path) as f:
        data = yaml.safe_load(f)

    project = data["project"]
    species_data = data["species"]

    species = {}

    for name, sp in species_data.items():
        species[name] = SpeciesConfig(
            name=name,
            genome=Path(sp["genome"]) if "genome" in sp else None,
            annotations=[Path(x) for x in sp.get("annotations", [])],
            proteins=[Path(x) for x in sp.get("proteins", [])],
        )

    return ProjectConfig(
        hal=Path(project["hal"]),
        species=species,
        references=project.get("references", []),
    )
