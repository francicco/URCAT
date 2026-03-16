from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.loci.comparative_loci import build_comparative_locus


def test_build_comparative_locus():
    hel_locus = SpeciesLocus(
        locus_id="hel_locus1",
        species="Heliconius",
        seqid="chr1",
        start=100,
        end=500,
        strand="+",
        transcripts=["tx1"]
    )

    dry_locus = SpeciesLocus(
        locus_id="dry_locus1",
        species="Dryas",
        seqid="chr1",
        start=120,
        end=480,
        strand="+",
        transcripts=["txA"]
    )

    projections = [
        ProjectedTranscript(
            species="Dryas",
            seqid="chr1",
            strand="+",
            source_species="Heliconius",
            source_transcript="tx1",
            exons=[(130, 470)],
            coverage=1.0,
        )
    ]

    species_loci = {
        "Heliconius": [hel_locus],
        "Dryas": [dry_locus],
    }

    clocus = build_comparative_locus(
        "clocus1",
        "Heliconius",
        "tx1",
        projections,
        species_loci,
    )

    assert clocus.primary["Dryas"] == "dry_locus1"
    assert clocus.species_count == 1
