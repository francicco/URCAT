from comparative_annotator.models.transcript import CandidateTranscript


def parse_attributes(attr):

    out = {}

    for item in attr.split(";"):
        if "=" not in item:
            continue
        k, v = item.split("=", 1)
        out[k] = v

    return out


def load_gff3(path, species):

    transcripts = {}

    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            seqid, _, ftype, start, end, _, strand, _, attr = line.strip().split("\t")

            attrs = parse_attributes(attr)

            start = int(start)
            end = int(end)

            if ftype.lower() in ["mrna", "transcript"]:

                tid = attrs["ID"]

                transcripts[tid] = CandidateTranscript(
                    transcript_id=tid,
                    species=species,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                )

            elif ftype.lower() == "exon":

                parent = attrs["Parent"]

                if parent in transcripts:
                    transcripts[parent].exons.append((start, end))

            elif ftype.lower() == "cds":

                parent = attrs["Parent"]

                if parent in transcripts:
                    transcripts[parent].cds.append((start, end))

    for tx in transcripts.values():
        tx.finalize()

    return transcripts
