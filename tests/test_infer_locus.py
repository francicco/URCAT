from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.io.hal import HALAdapter, HALCommandResult
from comparative_annotator.pipeline.infer_locus import infer_comparative_locus


def fake_runner(cmd):
    out_bed = cmd[-1]

    # simple mocked projection: every exon goes to chr5 with same coordinates
    with open(out_bed, "w") as f:
        f.write("chr5\t999\t1050\ttx1\t0\t+\n")

    return HALCommandResult(0, "", "")


def test_infer_comparative_locus():
    seed = CandidateTranscript(
        transcript_id="tx1",
        species="Heliconius",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        source="test",
        exons=[(100, 150), (300, 350)],
        cds=[],
    )
    seed.finalize()

    dry_locus = SpeciesLocus(
        locus_id="dry_locus1",
        species="Dryas",
        seqid="chr5",
        start=900,
        end=1200,
        strand="+",
        transcripts=["txA"],
    )

    hal = HALAdapter("dummy.hal", runner=fake_runner)

    clocus = infer_comparative_locus(
        seed_transcript=seed,
        target_species="Dryas",
        hal_adapter=hal,
        species_loci={"Dryas": [dry_locus]},
    )

    assert "Dryas" in clocus.primary
    assert clocus.primary["Dryas"] == "dry_locus1"
