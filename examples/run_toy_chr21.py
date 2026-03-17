from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus

from comparative_annotator.projection.matching import nearest_species_locus
from comparative_annotator.projection.reconstruct import reconstruct_projected_transcripts

from comparative_annotator.projection.reporting import rank_candidate_loci_with_transcripts

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

    transcripts_by_species = {
        "Hmel": hmel,
        "Eisa": eisa,
        "Diul": diul,
    }

    seed = hmel["transcript:Hmel202001oG3.1"]

    hal = HALAdapter("data/3SpChr21.hal")

    for target in ["Eisa", "Diul"]:
        projected_exon_blocks = []

        for exon_start, exon_end in seed.exons:
            intervals = hal.project_interval(
                source_species=seed.species,
                target_species=target,
                seqid=seed.seqid,
                start=exon_start,
                end=exon_end,
                strand=seed.strand,
                source_transcript=seed.transcript_id,
            )
            print(f"Projected exon {exon_start}-{exon_end} to {target}:")
            print(intervals)
            projected_exon_blocks.append(intervals)

        projected_transcripts = reconstruct_projected_transcripts(
            seed,
            projected_exon_blocks,
        )

        print(f"Reconstructed projected transcripts for {target}:")
        print(projected_transcripts)

        for pt in projected_transcripts:
            nearest, dist = nearest_species_locus(pt, species_loci[target])
            print("projected transcript:", pt)
            print("nearest locus:", nearest)
            print("distance:", dist)

        clocus = infer_comparative_locus(
            seed_transcript=seed,
            target_species=target,
            hal_adapter=hal,
            species_loci=species_loci,
            transcripts_by_species=transcripts_by_species,
        )

        print(f"=== {seed.transcript_id} -> {target} ===")
        print("primary:", clocus.primary)
        print("alternatives:", clocus.alternatives)
        print("strand_conflicts:", clocus.strand_conflicts)
        print("missing:", clocus.missing_annotations)
        print("primary_transcripts:", clocus.primary_transcripts)
        print("alternative_transcripts:", clocus.alternative_transcripts)

        nearest_same, dist_same = nearest_species_locus(pt, species_loci[target], same_strand_only=True)
        nearest_any, dist_any = nearest_species_locus(pt, species_loci[target], same_strand_only=False)

        print("nearest same-strand locus:", nearest_same)
        print("same-strand distance:", dist_same)
        print("nearest any-strand locus:", nearest_any)
        print("any-strand distance:", dist_any)
        print("=========================================\n")


if __name__ == "__main__":
    main()
