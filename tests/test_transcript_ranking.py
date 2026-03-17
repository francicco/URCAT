from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.projection.transcript_ranking import choose_best_transcript_within_locus


def test_choose_best_transcript_within_locus():
    source = CandidateTranscript(
        transcript_id="src_tx",
        species="Hmel",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        source="test",
        exons=[(100, 150), (200, 250), (300, 350)],
        cds=[],
    )
    source.finalize()

    projected = ProjectedTranscript(
        species="Diul",
        seqid="Diul2100",
        strand="-",
        source_species="Hmel",
        source_transcript="src_tx",
        exons=[(1000, 1050), (1100, 1150), (1200, 1250)],
        coverage=3,
        source_exon_indices=[0, 1, 2],
        chain_orientation="forward",
    )

    tx1 = CandidateTranscript(
        transcript_id="tx1",
        species="Diul",
        seqid="Diul2100",
        start=950,
        end=1300,
        strand="-",
        source="test",
        exons=[(1000, 1050), (1100, 1150), (1200, 1250)],
        cds=[],
    )
    tx1.finalize()

    tx2 = CandidateTranscript(
        transcript_id="tx2",
        species="Diul",
        seqid="Diul2100",
        start=950,
        end=1300,
        strand="+",
        source="test",
        exons=[(1000, 1030)],
        cds=[],
    )
    tx2.finalize()

    locus = SpeciesLocus(
        locus_id="diul_locus1",
        species="Diul",
        seqid="Diul2100",
        start=900,
        end=1300,
        strand="-",
        transcripts=["tx1", "tx2"],
    )

    best, alternatives = choose_best_transcript_within_locus(
        projected=projected,
        source_transcript=source,
        locus=locus,
        target_transcripts_by_id={"tx1": tx1, "tx2": tx2},
    )

    assert best.transcript_id == "tx1"
