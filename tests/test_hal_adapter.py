from pathlib import Path
from comparative_annotator.io.hal import HALAdapter, HALCommandResult


def fake_runner_success_factory(expected_output: str):
    def _runner(cmd):
        out_bed = Path(cmd[-1])
        out_bed.write_text(expected_output)
        return HALCommandResult(
            returncode=0,
            stdout="",
            stderr="",
        )
    return _runner


def test_project_interval_parses_bed_output():
    fake_bed = "chr5\t99\t200\ttx1\t0\t+\n"

    adapter = HALAdapter(
        "dummy.hal",
        runner=fake_runner_success_factory(fake_bed),
    )

    projections = adapter.project_interval(
        source_species="Heliconius",
        target_species="Dryas",
        seqid="chr1",
        start=100,
        end=200,
        strand="+",
        source_transcript="tx1",
    )

    assert len(projections) == 1
    p = projections[0]

    assert p.species == "Dryas"
    assert p.seqid == "chr5"
    assert p.start == 100
    assert p.end == 200
    assert p.strand == "+"


def test_to_bed_line_converts_coordinates():
    line = HALAdapter._to_bed_line(
        seqid="chr1",
        start=101,
        end=250,
        name="tx1",
        strand="-",
    )

    assert line == "chr1\t100\t250\ttx1\t0\t-\n"
