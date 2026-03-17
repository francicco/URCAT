from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.projection.adjudication import choose_best_locus


def test_choose_best_locus():
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

    good = SpeciesLocus(
        locus_id="good",
        species="Dryas",
        seqid="chr5",
        start=900,
        end=1200,
        strand="+",
        transcripts=["txA"],
    )

    bad = SpeciesLocus(
        locus_id="bad",
        species="Dryas",
        seqid="chr5",
        start=1130,
        end=1140,
        strand="+",
        transcripts=["txB"],
    )

    best, alternatives = choose_best_locus(projected, source, [good, bad])

    assert best.locus_id == "good"
    assert len(alternatives) == 1
