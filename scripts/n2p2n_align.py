# Takes multi-FASTA input, finds largest ORF for each, translates each ORF to the corresponding amino acid sequence, 
# performs multiple alignment with MUSCLE, and then converts the gapped amino acid sequence back to the nucleotide
# sequence, preserving the gaps and the original codons 
# USAGE: python n2p2n.py [multi_FASTA file]

# Import libraries
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import shlex
import subprocess

# Function that finds all ORFs and their locations (from Biopython website tutorial)
def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer

# Find largest ORF for each sequence in input FASTA file (orf_list format is [start, end, strand, sequence])
table = 12
min_pro_len = 100
orfs = []
ids = []
seqs = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
	ids.append(record.description)
	seqs.append(record.seq)
	orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
	orf_lens = []
	for start,end,strand,pro in orf_list:
		orf_lens.append(end-start)
	orfs.append(orf_list[orf_lens.index(max(orf_lens))])

# Write ORFs to multiFASTA file
with open("orfs.fasta", "w") as handle:
	for orf,name in zip(orfs,ids):
		handle.write(">"+name+"\n"+orf[3]+"\n")

# Call MUSCLE protein aligner and align ORFs 
subprocess.call(shlex.split('/home/nickeener/tools/muscle -in orfs.fasta -out aligned_orfs.fasta'))

# Delete orfs.fasta file
subprocess.call(['rm', 'orfs.fasta'])

# Read MUSCLE FASTA output and restore original order
align_seqs = {}
for record in SeqIO.parse("aligned_orfs.fasta", "fasta"):
	align_seqs[record.description] = record.seq
sorted_seqs = []
for name in ids:
	for key in align_seqs.keys():
		if name == key:
			sorted_seqs.append(align_seqs[key])
		else:
			pass

# Delete alignment file
subprocess.call(['rm', 'aligned_orfs.fasta'])

# Convert back to nucleotide sequence while preserving gaps and original codon sequence
new_seqs = []
nuc_seqs = []
for prot_seq,nuc_seq, x in zip(sorted_seqs,seqs,range(len(orfs))):
	start = orfs[x][0]
	end = orfs[x][1]
	nuc_seq = nuc_seq[start:end]
	new_seq = ''
	count = 0
	new_count = 0
	for i in str(prot_seq):
		if i == '-':
			count += 1
		else:
			pass
	for i,k in zip(range(len(prot_seq)),range(0, len(nuc_seq)+3*count, 3)):
		k -= new_count*3
		if prot_seq[i] != '-':
			new_seq += nuc_seq[k:k+3]
		else:
			new_seq += '---'
			new_count += 1
	new_seq += nuc_seq[-3:]
	new_seqs.append(str(new_seq))

# Write new sequences to file
with open("ultraAligned.fasta", "w") as handle:
	for seq,name in zip(new_seqs,ids):
		handle.write(">"+name+"\n"+seq+"\n")




