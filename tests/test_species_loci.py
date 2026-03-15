from pathlib import Path

from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.loci.species_loci import build_species_loci


def test_build_species_loci_overlapping_transcripts(tmp_path: Path):
    gff = tmp_path / "overlap.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\texon\t300\t400\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\tmRNA\t350\t600\t.\t+\t.\tID=tx2;Parent=gene2",
                "chr1\tsrc\texon\t350\t450\t.\t+\t.\tParent=tx2",
                "chr1\tsrc\texon\t500\t600\t.\t+\t.\tParent=tx2",
            ]
        )
    )

    txs = list(load_gff3(gff, species="sp1").values())
    loci = build_species_loci(txs, species="sp1")

    assert len(loci) == 1
    assert loci[0].transcript_count == 2
    assert set(loci[0].transcripts) == {"tx1", "tx2"}


def test_build_species_loci_non_overlapping_transcripts(tmp_path: Path):
    gff = tmp_path / "nonoverlap.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\tmRNA\t500\t600\t.\t+\t.\tID=tx2;Parent=gene2",
                "chr1\tsrc\texon\t500\t600\t.\t+\t.\tParent=tx2",
            ]
        )
    )

    txs = list(load_gff3(gff, species="sp1").values())
    loci = build_species_loci(txs, species="sp1")

    assert len(loci) == 2
    assert loci[0].transcripts == ["tx1"]
    assert loci[1].transcripts == ["tx2"]


def test_build_species_loci_separate_strands(tmp_path: Path):
    gff = tmp_path / "strand.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t300\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t300\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\tmRNA\t150\t280\t.\t-\t.\tID=tx2;Parent=gene2",
                "chr1\tsrc\texon\t150\t280\t.\t-\t.\tParent=tx2",
            ]
        )
    )

    txs = list(load_gff3(gff, species="sp1").values())
    loci = build_species_loci(txs, species="sp1")

    assert len(loci) == 2


def test_build_species_loci_separate_seqids(tmp_path: Path):
    gff = tmp_path / "seqid.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t300\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t300\t.\t+\t.\tParent=tx1",
                "chr2\tsrc\tmRNA\t150\t280\t.\t+\t.\tID=tx2;Parent=gene2",
                "chr2\tsrc\texon\t150\t280\t.\t+\t.\tParent=tx2",
            ]
        )
    )

    txs = list(load_gff3(gff, species="sp1").values())
    loci = build_species_loci(txs, species="sp1")

    assert len(loci) == 2


def test_locus_id_is_assigned_to_transcripts(tmp_path: Path):
    gff = tmp_path / "assign.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t300\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t300\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\tmRNA\t200\t400\t.\t+\t.\tID=tx2;Parent=gene2",
                "chr1\tsrc\texon\t200\t400\t.\t+\t.\tParent=tx2",
            ]
        )
    )

    txs = list(load_gff3(gff, species="sp1").values())
    loci = build_species_loci(txs, species="sp1")

    assert len(loci) == 1
    for tx in txs:
        assert tx.attributes["locus_id"] == loci[0].locus_id
