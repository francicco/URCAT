from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.projection.matching import find_overlapping_species_loci


def test_matching_single_hit():

    pt = ProjectedTranscript(
        species="Dryas",
        seqid="chr5",
        strand="+",
        source_species="Heliconius",
        source_transcript="tx1",
        exons=[(1000,1050),(1100,1150)]
    )

    locus = SpeciesLocus(
        locus_id="dry_locus1",
        species="Dryas",
        seqid="chr5",
        start=950,
        end=1200,
        strand="+",
        transcripts=["txA"]
    )

    hits = find_overlapping_species_loci(pt,[locus])

    assert len(hits) == 1
    assert hits[0].locus_id == "dry_locus1"
