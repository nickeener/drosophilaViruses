# script to cut a sequence into n length segments (overlapping and nonoverlapping) for use in BloomTree
import sys
from Bio import Entrez, SeqIO

global n
n = 8836 # size of segments

# Retrive sequence and sequence length from entrez using the accession number
def getseq(ID):
	Entrez.email = "nickeener@gmail.com"
	handle = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", retmode="xml")
	record = Entrez.read(handle, "genbank")
	handle.close()
	return (record[0]["GBSeq_sequence"], record[0]["GBSeq_length"])

# cut into non-overlapping segments of length n
def nooverlap(seq, length):
	segs = []
	for i in range(length/n):
		segs.append(("id"+str(i), seq[i*n:i*n+n]))
	return segs

# cut into overlapping segments of length n
def overlap(seq, length):
	segs = []
	for i in range(length-n+1):
		segs.append(("id"+str(i), seq[i:i+n]))
	return segs

# fetch sequence using accession number, cut sequence in n length segments, and print each segment with its
# ID to stdout (use ">" to write to file)
data = getseq(sys.argv[2])
lowerseq = data[0]
length = int(data[1])
upperseq = ""
for i in lowerseq:
	upperseq += i.upper()
if sys.argv[1] == "-n":
	out = nooverlap(upperseq, length)
	for i in range(length/n):
		print (out[i][0]+"\t"+out[i][1])
elif sys.argv[1] == "-o":
	out = overlap(upperseq, length)
	for i in range(length-n+1):
		print (out[i][0]+"\t"+out[i][1])
else:
	print ("Usage: cut.py [options] <Genbank Accession Number>")
	print ("Options:\n\t-n = no overlaps\n\t-o = overlaps")

'''Drosophila virus accession numbers:
Twyford Virus: KP714075
Kallithea Virus: NC_033829
DXC: NC_004177 (Segment A), NC_004169 (Segment B)
La Jolla Virus: NC_027128'''