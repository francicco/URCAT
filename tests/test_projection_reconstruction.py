from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.models.projection import ProjectionInterval
from comparative_annotator.projection.reconstruct import reconstruct_projected_transcripts


def test_reconstruct_single_clean_chain():
    tx = CandidateTranscript(
        transcript_id="tx1",
        species="Heliconius",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        source="test",
        exons=[(100, 150), (200, 250), (300, 350)],
        cds=[],
    )
    tx.finalize()

    projected_exon_blocks = [
        [ProjectionInterval("Dryas", "chr5", 1000, 1050, "+", "Heliconius", "tx1")],
        [ProjectionInterval("Dryas", "chr5", 1100, 1150, "+", "Heliconius", "tx1")],
        [ProjectionInterval("Dryas", "chr5", 1200, 1250, "+", "Heliconius", "tx1")],
    ]

    chains = reconstruct_projected_transcripts(tx, projected_exon_blocks)

    assert len(chains) == 1
    assert chains[0].seqid == "chr5"
    assert chains[0].strand == "+"
    assert chains[0].exons == [(1000, 1050), (1100, 1150), (1200, 1250)]
    assert chains[0].exon_count == 3


def test_reconstruct_splits_when_source_exon_order_breaks():
    tx = CandidateTranscript(
        transcript_id="tx1",
        species="Heliconius",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        source="test",
        exons=[(100, 150), (200, 250), (300, 350)],
        cds=[],
    )
    tx.finalize()

    # exon 2 projects before exon 1 in target coordinates -> should split
    projected_exon_blocks = [
        [ProjectionInterval("Dryas", "chr5", 2000, 2050, "+", "Heliconius", "tx1")],
        [ProjectionInterval("Dryas", "chr5", 1000, 1050, "+", "Heliconius", "tx1")],
        [ProjectionInterval("Dryas", "chr5", 3000, 3050, "+", "Heliconius", "tx1")],
    ]

    chains = reconstruct_projected_transcripts(tx, projected_exon_blocks)

    assert len(chains) >= 2


def test_reconstruct_splits_different_targets():
    tx = CandidateTranscript(
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
    tx.finalize()

    projected_exon_blocks = [
        [ProjectionInterval("Dryas", "chr5", 1000, 1050, "+", "Heliconius", "tx1")],
        [ProjectionInterval("Dryas", "chr8", 2000, 2050, "+", "Heliconius", "tx1")],
    ]

    chains = reconstruct_projected_transcripts(tx, projected_exon_blocks)

    assert len(chains) == 2
    seqids = sorted([c.seqid for c in chains])
    assert seqids == ["chr5", "chr8"]


def test_reconstruct_handles_multiple_blocks_per_source_exon():
    tx = CandidateTranscript(
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
    tx.finalize()

    projected_exon_blocks = [
        [
            ProjectionInterval("Dryas", "chr5", 1000, 1050, "+", "Heliconius", "tx1"),
            ProjectionInterval("Dryas", "chr9", 5000, 5050, "+", "Heliconius", "tx1"),
        ],
        [
            ProjectionInterval("Dryas", "chr5", 1100, 1150, "+", "Heliconius", "tx1"),
            ProjectionInterval("Dryas", "chr9", 5100, 5150, "+", "Heliconius", "tx1"),
        ],
    ]

    chains = reconstruct_projected_transcripts(tx, projected_exon_blocks)

    assert len(chains) == 2
    seqids = sorted(c.seqid for c in chains)
    assert seqids == ["chr5", "chr9"]
