from Bio import SeqIO
import sys

for record in SeqIO.parse(sys.argv[1], "fasta"):
	print(record.id+"\t"+record.seq+"\n")

