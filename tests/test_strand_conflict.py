from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.loci.comparative_builder import build_comparative_locus_from_projection


def test_strand_conflict_candidate():
    projected = ProjectedTranscript(
        species="Diul",
        seqid="Diul2100",
        strand="-",
        source_species="Hmel",
        source_transcript="tx1",
        exons=[(1000, 1100)],
        coverage=1.0,
    )

    locus = SpeciesLocus(
        locus_id="diul_locus1",
        species="Diul",
        seqid="Diul2100",
        start=900,
        end=1200,
        strand="+",
        transcripts=["txA"],
    )

    clocus = build_comparative_locus_from_projection(
        seed_species="Hmel",
        seed_transcript="tx1",
        projected_transcripts=[projected],
        species_loci={"Diul": [locus]},
    )

    assert "Diul" in clocus.strand_conflicts
    assert clocus.strand_conflicts["Diul"] == ["diul_locus1"]
    assert "Diul" not in clocus.missing_annotations
