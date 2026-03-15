# URCAT
A comparative annotation algorithm to annotate clade level group of genomes

```
comparative-annotator/
├── pyproject.toml
├── README.md
├── .gitignore
├── src/
│   └── comparative_annotator/
│       ├── __init__.py
│       ├── cli.py
│       ├── config/
│       │   ├── __init__.py
│       │   ├── manifest.py
│       │   └── validation.py
│       ├── io/
│       │   ├── __init__.py
│       │   ├── gff3.py
│       │   ├── gtf.py
│       │   ├── hal.py
│       │   └── writers.py
│       ├── models/
│       │   ├── __init__.py
│       │   ├── transcript.py
│       │   ├── locus.py
│       │   └── consensus.py
│       ├── normalize/
│       │   ├── __init__.py
│       │   └── annotations.py
│       ├── loci/
│       │   ├── __init__.py
│       │   ├── species_loci.py
│       │   └── comparative_loci.py
│       ├── scoring/
│       │   ├── __init__.py
│       │   ├── local.py
│       │   ├── projection.py
│       │   ├── comparative.py
│       │   └── priors.py
│       ├── adjudication/
│       │   ├── __init__.py
│       │   └── iterative.py
│       ├── qc/
│       │   ├── __init__.py
│       │   ├── outliers.py
│       │   └── reports.py
│       └── utils/
│           ├── __init__.py
│           ├── intervals.py
│           ├── logging.py
│           └── stats.py
├── tests/
│   ├── test_manifest.py
│   ├── test_gff3_parser.py
│   ├── test_species_loci.py
│   └── data/
│       ├── tiny_manifest.yaml
│       └── tiny.gff3
├── examples/
│   └── toy_project/
│       └── manifest.yaml
└── workflow/
    └── Snakefile
```

What each area is for

`config/`
reads and validates the YAML manifest.

`io/`
parses external files and writes outputs.
This is where GFF3/GTF/HAL adapters live.

`models/`
defines your main classes:
	•	transcript
	•	species locus
	•	comparative locus
	•	consensus

`normalize/`
converts parsed annotations into the internal format.

`loci/`
builds local loci and comparative loci.

`scoring/`
contains the scoring functions, split into local, projection, comparative, and prior.

`adjudication/`
contains the iterative algorithm that picks the best candidates.

`qc/`
flags weird annotations and suspicious loci.

`utils/`
small helper code only.


### Minimal pyproject.toml

```
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "comparative-annotator"
version = "0.1.0"
description = "Comparative annotation adjudication framework"
readme = "README.md"
requires-python = ">=3.11"
authors = [
  {name = "Francesco Cicconardi"}
]
dependencies = [
  "pyyaml>=6.0",
  "pytest>=8.0",
]

[project.scripts]
comparative-annotator = "comparative_annotator.cli:main"

[tool.setuptools]
package-dir = {"" = "src"}

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

# comparative-annotator

A comparative framework for adjudicating gene models across related genomes using multiple annotation sources and whole-genome alignment.

## Development install

```bash
pip install -e .
```

pytest

## First files to create

You do not need to fill everything at once.

I would start with only these:

- `pyproject.toml`
- `src/comparative_annotator/__init__.py`
- `src/comparative_annotator/cli.py`
- `src/comparative_annotator/config/manifest.py`
- `src/comparative_annotator/models/transcript.py`
- `src/comparative_annotator/models/locus.py`
- `src/comparative_annotator/io/gff3.py`
- `src/comparative_annotator/loci/species_loci.py`
- `tests/test_manifest.py`
- `tests/test_gff3_parser.py`
- `tests/test_species_loci.py`

That is enough for the first milestone.

## Minimal CLI

For now, `cli.py` can just be:

```python
def main():
    print("comparative-annotator")
```
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

```bash
comparative-annotator examples/manifest.yaml
```

First real development milestones

Next steps should be:
	1.	unit tests for the GFF3 parser
	
	2.	unit tests for species locus builder
	
	3.	HAL adapter interface
	
	4.	comparative locus builder
	
	5.	iterative adjudicator









