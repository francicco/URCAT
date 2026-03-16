from comparative_annotator.io.gff3 import load_gff3
from comparative_annotator.io.hal import HALAdapter, HALError


def main():
    hmel = load_gff3("data/Hmel202001o.test.gff3", species="Hmel")
    hal = HALAdapter("data/3SpChr21.hal")

    for tx_id, tx in hmel.items():
        print(f"\n=== {tx_id} ===")
        for target in ["Eisa", "Diul"]:
            exon_results = []
            for exon_start, exon_end in tx.exons:
                try:
                    intervals = hal.project_interval(
                        source_species=tx.species,
                        target_species=target,
                        seqid=tx.seqid,
                        start=exon_start,
                        end=exon_end,
                        strand=tx.strand,
                        source_transcript=tx.transcript_id,
                    )
                    exon_results.append(len(intervals))
                except HALError:
                    exon_results.append(0)

            print(f"{target}: exon block counts = {exon_results}")


if __name__ == "__main__":
    main()
