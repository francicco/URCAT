from __future__ import annotations

from pathlib import Path

from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.models.transcript import CandidateTranscript


def load_transcripts_for_species(
    annotation_paths: dict[str, str],
    species: str,
) -> dict[str, CandidateTranscript]:
    """
    Load transcripts for one species from the explicit per-species annotation map.

    Missing species are allowed and return {}.
    """
    gff_path = annotation_paths.get(species)

    if not gff_path:
        return {}

    p = Path(gff_path)
    if not p.exists():
        raise FileNotFoundError(f"Annotation for species '{species}' not found: {p}")

    return load_gff3(str(p), species=species)


def load_all_transcripts(
    annotation_paths: dict[str, str],
    species_list: list[str],
) -> dict[str, dict[str, CandidateTranscript]]:
    """
    Load all available transcript annotations.
    Species without annotation return {}.
    """
    return {
        species: load_transcripts_for_species(annotation_paths, species)
        for species in species_list
    }


def build_all_species_loci(
    transcripts_by_species: dict[str, dict[str, CandidateTranscript]],
):
    """
    Build native loci per species from loaded transcripts.
    Species without annotation produce an empty locus list.
    """
    return {
        species: build_species_loci(list(txdict.values()), species=species)
        if txdict
        else []
        for species, txdict in transcripts_by_species.items()
    }


def species_has_annotation(
    annotation_paths: dict[str, str],
    species: str,
) -> bool:
    return species in annotation_paths and bool(annotation_paths[species])


def available_annotation_species(
    annotation_paths: dict[str, str],
    species_list: list[str],
) -> list[str]:
    return [species for species in species_list if species_has_annotation(annotation_paths, species)]
