from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.missing.consensus import (
    cluster_projected_transcripts,
    choose_missing_locus_strand,
    build_consensus_missing_transcript,
)


def test_choose_missing_locus_strand():
    pts = [
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="+",
            source_species="Hmel",
            source_transcript="tx1",
            exons=[(100, 150), (200, 250)],
            coverage=2,
            source_exon_indices=[0, 1],
            chain_score=2.0,
        ),
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="+",
            source_species="Diul",
            source_transcript="tx2",
            exons=[(102, 152), (198, 248)],
            coverage=2,
            source_exon_indices=[0, 1],
            chain_score=1.8,
        ),
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="-",
            source_species="Other",
            source_transcript="tx3",
            exons=[(105, 155), (205, 255)],
            coverage=2,
            source_exon_indices=[0, 1],
            chain_score=1.0,
        ),
    ]

    best_strand, strand_scores = choose_missing_locus_strand(pts)

    assert best_strand == "+"
    assert "+" in strand_scores
    assert "-" in strand_scores


def test_build_consensus_missing_transcript():
    pts = [
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="+",
            source_species="Hmel",
            source_transcript="tx1",
            exons=[(100, 150), (200, 250)],
            coverage=2,
            source_exon_indices=[0, 1],
            chain_score=2.0,
        ),
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="+",
            source_species="Diul",
            source_transcript="tx2",
            exons=[(102, 152), (198, 248)],
            coverage=2,
            source_exon_indices=[0, 1],
            chain_score=1.8,
        ),
    ]

    consensus = build_consensus_missing_transcript(pts, chosen_strand="+")

    assert consensus.species == "Eisa"
    assert consensus.seqid == "Eisa2100"
    assert consensus.strand == "+"
    assert len(consensus.exons) == 2
    assert consensus.exons[0][0] in (101, 102)
    assert consensus.exons[0][1] in (151, 152)


def test_cluster_projected_transcripts():
    pts = [
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="+",
            source_species="Hmel",
            source_transcript="tx1",
            exons=[(100, 150)],
        ),
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="-",
            source_species="Diul",
            source_transcript="tx2",
            exons=[(120, 180)],
        ),
        ProjectedTranscript(
            species="Eisa",
            seqid="Eisa2100",
            strand="+",
            source_species="Other",
            source_transcript="tx3",
            exons=[(1000, 1100)],
        ),
    ]

    clusters = cluster_projected_transcripts(pts, max_gap=0)

    assert len(clusters) == 2
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 1
