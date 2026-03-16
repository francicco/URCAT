from comparative_annotator.models.transcript import CandidateTranscript
from comparative_annotator.models.locus import SpeciesLocus
from comparative_annotator.io.hal import HALAdapter, HALCommandResult
from comparative_annotator.projection.reciprocal import validate_reciprocal_projection


def fake_runner(cmd):
    out_bed = cmd[-1]

    # Always project to chr5:1000-1050 (+)
    with open(out_bed, "w") as f:
        f.write("chr5\t999\t1050\ttx1\t0\t+\n")

    return HALCommandResult(0, "", "")


def test_validate_reciprocal_projection():
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

    dry_rep = CandidateTranscript(
        transcript_id="dry_tx1",
        species="Dryas",
        seqid="chr5",
        start=900,
        end=1200,
        strand="+",
        source="test",
        exons=[(1000, 1050), (1100, 1150)],
        cds=[],
    )
    dry_rep.finalize()

    hel_locus = SpeciesLocus(
        locus_id="hel_locus1",
        species="Heliconius",
        seqid="chr1",
        start=50,
        end=500,
        strand="+",
        transcripts=["tx1"],
    )

    dry_locus = SpeciesLocus(
        locus_id="dry_locus1",
        species="Dryas",
        seqid="chr5",
        start=900,
        end=1200,
        strand="+",
        transcripts=["dry_tx1"],
    )

    species_loci = {
        "Heliconius": [hel_locus],
        "Dryas": [dry_locus],
    }

    transcript_lookup = {
        ("Dryas", "dry_locus1"): dry_rep,
    }

    def representative_transcript_getter(species, locus_id):
        return transcript_lookup[(species, locus_id)]

    hal = HALAdapter("dummy.hal", runner=fake_runner)

    result = validate_reciprocal_projection(
        seed_transcript=seed,
        source_species_locus_id="hel_locus1",
        target_species="Dryas",
        hal_adapter=hal,
        species_loci=species_loci,
        representative_transcript_getter=representative_transcript_getter,
    )

    assert result.target_species == "Dryas"
    assert result.forward_primary_locus == "dry_locus1"
    assert result.classification in {
        "reciprocal_ortholog",
        "nonreciprocal_candidate",
        "missing_reverse_match",
    }
