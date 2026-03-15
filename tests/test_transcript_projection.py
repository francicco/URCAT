from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.io.hal import HALAdapter, HALCommandResult


def fake_runner(cmd):

    out_bed = cmd[-1]

    with open(out_bed, "w") as f:
        f.write("chr5\t99\t150\ttx1\t0\t+\n")

    return HALCommandResult(0, "", "")


def test_project_transcript():

    tx = CandidateTranscript(
        transcript_id="tx1",
        species="Heliconius",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        source="test",
        exons=[(100,150),(300,350)],
        cds=[]
    )

    adapter = HALAdapter("dummy.hal", runner=fake_runner)

    projected = adapter.project_transcript(tx,"Dryas")

    assert len(projected) == 1
    p = projected[0]

    assert p.species == "Dryas"
    assert p.seqid == "chr5"
    assert p.exon_count == 2
