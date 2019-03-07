import sys
import subprocess
from Bio import SeqIO

records = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
	records.append([str(record.id), str(record.seq)])
print(records[0][1])
#for seq
#subprocess.call(['search-fastq.sh', ])