from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--round-id", type=int, required=True)
    parser.add_argument("--reference-species", required=True)
    parser.add_argument("--target-species", required=True)
    parser.add_argument("--batch-id", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out = {
        "round_id": args.round_id,
        "reference_species": args.reference_species,
        "target_species": args.target_species,
        "batch_id": args.batch_id,
        "status": "done",
        "results": [],
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
