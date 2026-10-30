import argparse 
import gzip


parser = argparse.ArgumentParser(description= "Map transcripts to genes")
parser.add_argument("--raw", required=True, help="GENCODE transcripts FASTA")
parser.add_argument("--transcriptome", required=True, help="output FASTA with headers consisting only of transcript ids")
parser.add_argument("--tx2gene", required=True, help="output TSV consisting of mapping between transcript ids and gene symbols")

args = parser.parse_args()

with gzip.open(args.raw, "rt") as raw, gzip.open(args.transcriptome, "wt") as transcriptome, open(args.tx2gene, "wt") as tx2gene:
    for line in raw:
        if line.startswith(">"):
            col_names = line[1:].split("|")
            transcript_id = col_names[0]
            gene_symbol = col_names[5]
            transcriptome.write(f">{transcript_id}\n")
            tx2gene.write(f"{transcript_id}\t{gene_symbol}\n")
        else:
            transcriptome.write(line)

print("Reformatted transcriptome!")