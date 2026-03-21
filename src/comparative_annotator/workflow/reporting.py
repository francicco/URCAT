from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable


def read_json(path: str | Path):
    path = Path(path)
    with open(path) as fh:
        return json.load(fh)


def write_json(path: str | Path, obj) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def normalize_scalar(value):
    if isinstance(value, (list, tuple)):
        return ",".join("" if x is None else str(x) for x in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True)
    return value


def write_tsv(path: str | Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        with open(path, "w") as fh:
            if columns:
                fh.write("\t".join(columns) + "\n")
        return

    if columns is None:
        columns = sorted({k for row in rows for k in row.keys()})

    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            vals = []
            for col in columns:
                v = normalize_scalar(row.get(col))
                vals.append("" if v is None else str(v))
            fh.write("\t".join(vals) + "\n")
