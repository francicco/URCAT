from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.projection.scoring import score_projected_transcript_against_locus


def test_score_projected_transcript_against_locus():
    projected = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000, 1050), (1100, 1150)],
        coverage=1.0,
    )

    locus = SpeciesLocus(
        locus_id="dry_locus1",
        species="Dryas",
        seqid="chr5",
        start=950,
        end=1200,
        strand="+",
        transcripts=["txA"],
    )

    score = score_projected_transcript_against_locus(projected, locus)

    assert score.species == "Dryas"
    assert score.locus_id == "dry_locus1"
    assert score.total_score > 0
    assert 0 <= score.overlap_fraction <= 1


def test_score_is_lower_for_poor_overlap():
    projected = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000, 1050), (1100, 1150)],
        coverage=1.0,
    )

    good_locus = SpeciesLocus(
        locus_id="good",
        species="Dryas",
        seqid="chr5",
        start=950,
        end=1200,
        strand="+",
        transcripts=["txA"],
    )

    bad_locus = SpeciesLocus(
        locus_id="bad",
        species="Dryas",
        seqid="chr5",
        start=1140,
        end=1160,
        strand="+",
        transcripts=["txB"],
    )

    good_score = score_projected_transcript_against_locus(projected, good_locus)
    bad_score = score_projected_transcript_against_locus(projected, bad_locus)

    assert good_score.total_score > bad_score.total_score
