from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.projection.matching import (
    match_projected_transcript_to_loci,
    classify_unmatched_projection,
)


def test_match_single_locus():
    projected = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000, 1050), (1100, 1150)],
    )

    loci = [
        SpeciesLocus(
            locus_id="dry_locus1",
            species="Dryas",
            seqid="chr5",
            start=900,
            end=1200,
            strand="+",
            transcripts=["txA"],
        )
    ]

    matches = match_projected_transcript_to_loci(projected, loci)

    assert len(matches) == 1
    assert matches[0].classification == "ortholog_candidate"
    assert matches[0].locus_id == "dry_locus1"


def test_match_multiple_loci():
    projected = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000, 1050), (1100, 1150)],
    )

    loci = [
        SpeciesLocus(
            locus_id="dry_locus1",
            species="Dryas",
            seqid="chr5",
            start=950,
            end=1080,
            strand="+",
            transcripts=["txA"],
        ),
        SpeciesLocus(
            locus_id="dry_locus2",
            species="Dryas",
            seqid="chr5",
            start=1090,
            end=1300,
            strand="+",
            transcripts=["txB"],
        ),
    ]

    matches = match_projected_transcript_to_loci(projected, loci)

    assert len(matches) == 2
    assert all(m.classification == "paralog_candidate" for m in matches)


def test_unmatched_projection():
    projected = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000, 1050), (1100, 1150)],
    )

    loci = [
        SpeciesLocus(
            locus_id="dry_locus1",
            species="Dryas",
            seqid="chr8",
            start=900,
            end=1200,
            strand="+",
            transcripts=["txA"],
        )
    ]

    cls = classify_unmatched_projection(projected, loci)

    assert cls == "missing_annotation_candidate"
