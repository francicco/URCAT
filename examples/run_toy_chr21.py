from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci
from comparative_annotator.io.hal import HALAdapter
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus

from comparative_annotator.projection.matching import nearest_species_locus
from comparative_annotator.projection.reconstruct import reconstruct_projected_transcripts
from comparative_annotator.projection.reporting import rank_candidate_loci_with_transcripts

from comparative_annotator.missing.consensus import (
    cluster_projected_transcripts,
    choose_missing_locus_strand,
    build_consensus_missing_transcript,
)

from comparative_annotator.io.gff3_writer import write_full_merged_gff3


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

    # collect new URCAT loci for Eisa here
    eisa_consensus = []

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
            print("projected transcript:", pt)

            nearest, dist = nearest_species_locus(pt, species_loci[target])
            print("nearest locus:", nearest)
            print("distance:", dist)

            nearest_same, dist_same = nearest_species_locus(
                pt, species_loci[target], same_strand_only=True
            )
            nearest_any, dist_any = nearest_species_locus(
                pt, species_loci[target], same_strand_only=False
            )

            print("nearest same-strand locus:", nearest_same)
            print("same-strand distance:", dist_same)
            print("nearest any-strand locus:", nearest_any)
            print("any-strand distance:", dist_any)

            candidate_rows = rank_candidate_loci_with_transcripts(
                projected=pt,
                source_transcript=seed,
                species_loci=species_loci[target],
                target_transcripts_by_id=transcripts_by_species[target],
            )

            if candidate_rows:
                print("candidate loci ranking:")
                for i, row in enumerate(candidate_rows, start=1):
                    print(
                        f"  {i}. {row['locus_id']} "
                        f"relation={row['relation']} "
                        f"distance={row['distance']} "
                        f"locus_score={row['locus_score']:.3f} "
                        f"overlap={row['overlap_fraction']:.3f} "
                        f"proj_score={row['projection_score']:.3f} "
                        f"same_strand={row['same_strand']} "
                        f"best_tx={row['best_transcript_id']} "
                        f"best_tx_score={row['best_transcript_score']}"
                    )

        # Missing-locus consensus diagnostics + collect new Eisa loci
        if target == "Eisa" and projected_transcripts:
            print("missing-locus consensus diagnostics:")

            clusters = cluster_projected_transcripts(projected_transcripts, max_gap=0)
            print(f"number of clusters: {len(clusters)}")

            for i, cluster in enumerate(clusters, start=1):
                print(f"cluster {i}:")
                for pt in cluster:
                    print("  support transcript:", pt)

                best_strand, strand_scores = choose_missing_locus_strand(cluster)
                consensus = build_consensus_missing_transcript(cluster, best_strand)

                print("  chosen strand:", best_strand)
                print("  strand scores:")
                for strand, info in strand_scores.items():
                    print(
                        f"    {strand}: "
                        f"support_count={info['support_count']} "
                        f"total_chain_score={info['total_chain_score']:.3f} "
                        f"mean_exon_recovery={info['mean_exon_recovery']:.3f} "
                        f"score={info['score']:.3f}"
                    )

                print("  consensus transcript:", consensus)

                if consensus is not None:
                    eisa_consensus.append(consensus)

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
        print("=========================================\n")

    # write merged Eisa annotation at the end
    urcat_consensus_by_species = {
        "Eisa": eisa_consensus,
        "Hmel": [],
        "Diul": [],
    }

    write_full_merged_gff3(
        output_path="all_species_merged_urcat.gff3",
        original_transcripts_by_species=transcripts_by_species,
        urcat_consensus_by_species=urcat_consensus_by_species,
    )

    print("Wrote merged annotation to all_species_merged_urcat.gff3")


if __name__ == "__main__":
    main()
