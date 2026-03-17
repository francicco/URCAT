from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.loci.comparative_builder import build_comparative_locus_from_projection


def test_build_comparative_locus():
    source = CandidateTranscript(
        transcript_id="tx1",
        species="Heliconius",
        seqid="chr1",
        start=100,
        end=300,
        strand="+",
        source="test",
        exons=[(100, 150), (200, 250)],
        cds=[],
    )
    source.finalize()

    projected = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000, 1050), (1100, 1150)],
        coverage=1.0,
        source_exon_indices=[0, 1],
    )

    locus = SpeciesLocus(
        locus_id="dry_locus1",
        species="Dryas",
        seqid="chr5",
        start=900,
        end=1200,
        strand="+",
        transcripts=["txA"]
    )

    species_loci = {"Dryas": [locus]}

    clocus = build_comparative_locus_from_projection(
        "Heliconius",
        "tx1",
        [projected],
        species_loci,
        source_transcript=source,
        transcripts_by_species={},
    )

    assert "Dryas" in clocus.primary
    assert clocus.primary["Dryas"] == "dry_locus1"
