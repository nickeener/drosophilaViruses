# Creates n number of permuations of each seqeunce in a multi-FASTA file and saves the result (for each virus) in a separate
# FASTA file

# Import libraries
import sys
import subprocess
from Bio import SeqIO
import pdb

# Set number of permutations
n = 100
seqdir = '/home/nickeener/projects/drosophilaViruses/viralSeqs/'

# Read multi-FASTA containing 50 viral sequences and save description and sequence of each
names = []
seqs = []
for record in SeqIO.parse(seqdir+'drosophilaViruses.fasta', 'fasta'):
	names.append(record.description)
	seqs.append(record.seq)

# Create temporary directory to write individual sequences to FASTA files to and write the files
subprocess.call(['rm', '-r', seqdir+'tmp']) # Deletes tmp directory if it exists
subprocess.call(['mkdir', seqdir+'tmp'])
for name,seq in zip(names,seqs):
	with open(seqdir+'tmp/'+name.split()[0]+'.fasta', 'w') as file:
		file.write('>'+name+'\n'+str(seq))

# Call shuffleseq on each new FASTA file and create multiFASTA files containing n permutations of each virus
subprocess.call(['mkdir', seqdir+'permutations'])
for name in names:
	subprocess.call(['shuffleseq', '-sequence', seqdir+'tmp/'+name.split()[0]+'.fasta', '-outseq', seqdir+'permutations/'+name.split()[0]+'.fasta', '-shuffle', str(n)])

# Delete temporary files
subprocess.call(['rm', '-r', seqdir+'tmp'])
