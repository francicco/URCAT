from pathlib import Path
from tempfile import NamedTemporaryFile

from comparative_annotator.io.hal import HALAdapter, HALCommandResult

 def project_transcript(
     self,
     transcript,
     target_species: str,
 ) -> list[ProjectedTranscript]:
     """
     Project each exon of a transcript independently and group the results
     into projected transcript candidates.

     Version 1 grouping is deliberately simple:
     projected exons are grouped by (seqid, strand).
     """
     exon_projections = []

     for exon_start, exon_end in transcript.exons:
         intervals = self.project_interval(
             source_species=transcript.species,
             target_species=target_species,
             seqid=transcript.seqid,
             start=exon_start,
             end=exon_end,
             strand=transcript.strand,
             source_transcript=transcript.transcript_id,
         )
         exon_projections.append(intervals)

     flat: list = []
     for block in exon_projections:
         flat.extend(block)

     if not flat:
         return []

     grouped: dict[tuple[str, str], ProjectedTranscript] = {}

     for proj in flat:
         key = (proj.seqid, proj.strand)

         if key not in grouped:
             grouped[key] = ProjectedTranscript(
                 species=proj.species,
                 seqid=proj.seqid,
                 strand=proj.strand,
                 source_species=proj.source_species,
                 source_transcript=proj.source_transcript,
             )

         grouped[key].add_exon(proj.start, proj.end)

     return list(grouped.values())


def fake_runner_success_factory(expected_output: str):
    def _runner(cmd):
        # last argument is output BED path
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
    assert p.source_species == "Heliconius"
    assert p.source_transcript == "tx1"


def test_to_bed_line_converts_coordinates():
    line = HALAdapter._to_bed_line(
        seqid="chr1",
        start=101,
        end=250,
        name="tx1",
        strand="-",
    )

    assert line == "chr1\t100\t250\ttx1\t0\t-\n"
