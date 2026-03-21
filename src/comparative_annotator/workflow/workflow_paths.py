from pathlib import Path


def get_round_dir(workdir: str, round_id: int) -> Path:
    return Path(workdir) / "rounds" / f"round_{round_id:03d}"


def get_reference_dir(workdir: str, round_id: int, reference_species: str) -> Path:
    return get_round_dir(workdir, round_id) / f"ref_{reference_species}"


def get_target_dir(workdir: str, round_id: int, reference_species: str, target_species: str) -> Path:
    return get_reference_dir(workdir, round_id, reference_species) / f"target_{target_species}"


def get_summary_path(workdir: str, round_id: int) -> Path:
    return get_round_dir(workdir, round_id) / "summary.tsv"


def get_decision_path(workdir: str, round_id: int) -> Path:
    return get_round_dir(workdir, round_id) / "post_round_decision.json"


def get_edge_evidence_path(workdir: str, round_id: int, reference_species: str, target_species: str) -> Path:
    return get_target_dir(workdir, round_id, reference_species, target_species) / "edge_evidence.json"
