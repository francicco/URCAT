import argparse
from comparative_annotator.config.manifest import load_manifest

def main():

    parser = argparse.ArgumentParser(
        description="Comparative annotation adjudication framework"
    )

    parser.add_argument(
        "manifest",
        help="Project manifest YAML"
    )

    args = parser.parse_args()

    cfg = load_manifest(args.manifest)

    print("Loaded project:")
    print(cfg)

