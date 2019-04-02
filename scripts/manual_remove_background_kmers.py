
import sys
from Bio import SeqIO
import subprocess
import pandas as pd
import pdb

# Set study accession number and kmer value to a variable
study = sys.argv[1]
kmer = str(20)

# Set absolute path to query file to a variable
query = '/home/nickeener/projects/drosophilaViruses/viralSeqs/'+sys.argv[3]+'.fasta'

# Set drive type (0 for main drive, 1 for external drive)
drive = int(sys.argv[2])

# Read viral FASTA files and set sequence to variable
for record in SeqIO.parse('/home/nickeener/projects/drosophilaViruses/viralSeqs/DCV.fasta', "fasta"):
	dcv = str(record.seq)
for record in SeqIO.parse('/home/nickeener/projects/drosophilaViruses/viralSeqs/nora_virus.fasta', "fasta"):
	nora = str(record.seq)
for record in SeqIO.parse('/home/nickeener/projects/drosophilaViruses/viralSeqs/sigma_virus.fasta', "fasta"):
	sigma = str(record.seq)
for record in SeqIO.parse('/home/nickeener/projects/drosophilaViruses/viralSeqs/wolbachia.fasta', "fasta"):
	wolbachia = str(record.seq)
for record in SeqIO.parse(query, "fasta"):
	twyford = str(record.seq)

# Create list with every kmer in each viral sequence
dcv_kmers = []
nora_kmers = []
sigma_kmers = []
wolbachia_kmers = []
twyford_kmers = []
for i in range(len(dcv)):
	dcv_kmers.append(dcv[i:i+19])
for i in range(len(nora)):
	nora_kmers.append(nora[i:i+19])
for i in range(len(sigma)):
	sigma_kmers.append(sigma[i:i+19])
for i in range(len(wolbachia)):
	wolbachia_kmers.append(wolbachia[i:i+19])	
for i in range(len(twyford)):
	twyford_kmers.append(twyford[i:i+19])

# Remove duplicates
print(len(dcv_kmers))
dcv_kmers = list(dict.fromkeys(dcv_kmers))
print(len(dcv_kmers))
print(len(nora_kmers))
nora_kmers = list(dict.fromkeys(nora_kmers))
print(len(nora_kmers))
print(len(sigma_kmers))
sigma_kmers = list(dict.fromkeys(sigma_kmers))
print(len(sigma_kmers))
print(len(wolbachia_kmers))
wolbachia_kmers = list(dict.fromkeys(wolbachia_kmers))
print(len(wolbachia_kmers))
print(len(twyford_kmers))
twyford_kmers = list(dict.fromkeys(twyford_kmers))
print(len(twyford_kmers))

# Calculate kmers shared between Twyford and other viruses
print(len(list(set(twyford_kmers) & set(dcv_kmers))))
print(len(list(set(twyford_kmers) & set(nora_kmers))))
print(len(list(set(twyford_kmers) & set(wolbachia_kmers))))
print(len(list(set(twyford_kmers) & set(sigma_kmers))))


