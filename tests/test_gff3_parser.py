from pathlib import Path

from comparative_annotator.io.gff3 import load_gff3


def test_parse_simple_gff3(tmp_path: Path):
    gff = tmp_path / "tiny.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\texon\t300\t500\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\tCDS\t120\t180\t.\t+\t0\tParent=tx1",
                "chr1\tsrc\tCDS\t320\t450\t.\t+\t2\tParent=tx1",
            ]
        )
    )

    txs = load_gff3(gff, species="sp1")
    assert "tx1" in txs

    tx = txs["tx1"]
    assert tx.seqid == "chr1"
    assert tx.strand == "+"
    assert tx.exon_count == 2
    assert tx.intron_count == 1
    assert tx.exons == [(100, 200), (300, 500)]
    assert tx.cds == [(120, 180), (320, 450)]
    assert tx.intron_chain == [(201, 299)]


def test_exons_before_transcript_line(tmp_path: Path):
    gff = tmp_path / "unordered.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\texon\t300\t400\t.\t+\t.\tParent=tx1",
                "chr1\tsrc\tCDS\t120\t180\t.\t+\t0\tParent=tx1",
                "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene1",
            ]
        )
    )

    txs = load_gff3(gff, species="sp1")
    tx = txs["tx1"]

    assert tx.exon_count == 2
    assert tx.exons == [(100, 200), (300, 400)]
    assert tx.cds == [(120, 180)]


def test_multiple_parents(tmp_path: Path):
    gff = tmp_path / "multi_parent.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx2;Parent=gene1",
                "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=tx1,tx2",
                "chr1\tsrc\texon\t300\t500\t.\t+\t.\tParent=tx1,tx2",
            ]
        )
    )

    txs = load_gff3(gff, species="sp1")

    assert txs["tx1"].exons == [(100, 200), (300, 500)]
    assert txs["tx2"].exons == [(100, 200), (300, 500)]


def test_transcript_without_cds(tmp_path: Path):
    gff = tmp_path / "no_cds.gff3"
    gff.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t100\t200\t.\t-\t.\tID=tx1;Parent=gene1",
                "chr1\tsrc\texon\t100\t200\t.\t-\t.\tParent=tx1",
            ]
        )
    )

    txs = load_gff3(gff, species="sp1")
    tx = txs["tx1"]

    assert tx.exon_count == 1
    assert tx.cds_length == 0
    assert tx.cds == []
