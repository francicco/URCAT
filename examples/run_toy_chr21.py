from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus


def main():
    hmel = load_gff3("data/Hmel202001o.test.gff3", species="Hmel")
    eisa = load_gff3("data/Eisa2100.test.gff3", species="Eisa")
    diul = load_gff3("data/Diul2100.test.gff3", species="Diul")

    hmel_loci = build_species_loci(list(hmel.values()), species="Hmel")
    eisa_loci = build_species_loci(list(eisa.values()), species="Eisa")
    diul_loci = build_species_loci(list(diul.values()), species="Diul")

    species_loci = {
        "Hmel": hmel_loci,
        "Eisa": eisa_loci,
        "Diul": diul_loci,
    }

    seed = hmel["transcript:Hmel202001oG1.1"]

    hal = HALAdapter("data/3SpChr21.hal")

    for target in ["Eisa", "Diul"]:
        clocus = infer_comparative_locus(
            seed_transcript=seed,
            target_species=target,
            hal_adapter=hal,
            species_loci=species_loci,
        )

        print(f"=== {seed.transcript_id} -> {target} ===")
        print("primary:", clocus.primary)
        print("alternatives:", clocus.alternatives)
        print("missing:", clocus.missing_annotations)


if __name__ == "__main__":
    main()
