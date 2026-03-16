from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus


def main():
    hmel = load_gff3("/mnt/data/Hmel202001o.test.gff3", species="Hmel202001o")
    eisa = load_gff3("/mnt/data/Eisa2100.test.gff3", species="Eisa2100")
    diul = load_gff3("/mnt/data/Diul2100.test.gff3", species="Diul2100")

    hmel_loci = build_species_loci(list(hmel.values()), species="Hmel202001o")
    eisa_loci = build_species_loci(list(eisa.values()), species="Eisa2100")
    diul_loci = build_species_loci(list(diul.values()), species="Diul2100")

    species_loci = {
        "Hmel202001o": hmel_loci,
        "Eisa2100": eisa_loci,
        "Diul2100": diul_loci,
    }

    seed = hmel["transcript:Hmel202001oG1.1"]

    hal = HALAdapter("/mnt/data/3SpChr21.hal")

    for target in ["Eisa2100", "Diul2100"]:
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
