import sys
from Bio import Entrez, SeqIO


# Retrive sequence and sequence length from entrez using the accession number
def getseq(ID):
	Entrez.email = "nickeener@gmail.com"
	handle = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", retmode="xml")
	record = Entrez.read(handle, "genbank")
	handle.close()
	return (record[0]["GBSeq_sequence"])

seq = getseq(sys.argv[1])
print seq