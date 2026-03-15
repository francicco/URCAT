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
