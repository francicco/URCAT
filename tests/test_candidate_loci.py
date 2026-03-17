from comparative_annotator.models.projected_transcript import ProjectedTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.projection.matching import find_candidate_species_loci


def test_find_candidate_species_loci_includes_overlap_and_flanks():
    projected = ProjectedTranscript(
        species="Diul",
        seqid="Diul2100",
        strand="-",
        source_species="Hmel",
        source_transcript="tx1",
        exons=[(1000, 1100)],
    )

    loci = [
        SpeciesLocus("up2", "Diul", "Diul2100", 100, 200, "+", ["a"]),
        SpeciesLocus("up1", "Diul", "Diul2100", 700, 900, "+", ["b"]),
        SpeciesLocus("ov", "Diul", "Diul2100", 1050, 1300, "+", ["c"]),
        SpeciesLocus("down1", "Diul", "Diul2100", 1400, 1500, "+", ["d"]),
        SpeciesLocus("down2", "Diul", "Diul2100", 1600, 1700, "+", ["e"]),
        SpeciesLocus("down3", "Diul", "Diul2100", 2000, 2100, "+", ["f"]),
    ]

    candidates = find_candidate_species_loci(projected, loci, n_flank=2)
    ids = [x.locus_id for x in candidates]

    assert "ov" in ids
    assert "up1" in ids
    assert "up2" in ids
    assert "down1" in ids
    assert "down2" in ids
    assert "down3" not in ids
