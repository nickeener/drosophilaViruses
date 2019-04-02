# Creates file with number of unique kmers for every run in a set of reads
# USAGE: python kmer_counts.py [Study Accession Number] [Drive Type]

# Import libraries
import sys
import subprocess
import pandas as pd
import pdb

# Set study accession number and drive type to variables
study = sys.argv[1]
drive = int(sys.argv[2])

# Get list of all files/directories in the mapping directory of interest
if drive == 0:
	ls = subprocess.Popen(['ls', '/home/nickeener/projects/drosophilaViruses/mapping/'+study], stdout=subprocess.PIPE)
else:
	ls = subprocess.Popen(['ls', '/media/nickeener/External_Drive/'+study], stdout=subprocess.PIPE)
output = ls.stdout.read()

# Convert output string into list containing only the SRA run accession number of each file
reads = []
for i in range(len(output)):
	if output[i] == 'S' or output[i] == 'E' and output[i+1] == 'R':
		reads.append(output[i:i+21])

# Assign each mate to a different list
pairs1 = []
pairs2 = []
for i in range(len(reads)):
	if i%2 == 0:
		pairs1.append(reads[i])
	else:
		pairs2.append(reads[i])

# Get list of each run acession number
run_acc = []
for pair in pairs1:
	run_acc.append(pair[0:10])

# For each run, calculate the number of unique kmers that occur and store in file 
count_file = open("kmer_counts.txt", "w")
for pair1,pair2,run in zip(pairs1, pairs2, run_acc):
	print("jellyfish count "+run)
	# Create a file with the counts of every 20mer in the read data
	subprocess.call(['./kmer_counts.sh', study, run, pair1, pair2])
	# Extract data from kmer count file in chunks (too big)
	kmer_count = []
	chunk_number = 1
	for chunk in pd.read_csv("/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/dump.fa", names=['1'], chunksize=10000):
		data = chunk
		# Testing purposes
		print("Chunk Number: "+str(chunk_number))
		chunk_number += 1
		# Assign kmers to one list and counts to another (in order)
		chunk_kmers = []
		chunk_counts = []
		for i in range(len(data['1'])):
			if i%2 == 0:
				chunk_counts.append(data['1'][i])
			else: 
				chunk_kmers.append(data['1'][i])
		# Remove '>' from each count entry and convert to integers
		for i in range(len(chunk_counts)):
			chunk_counts[i] = int(chunk_counts[i][1:])
		# Remove all entries with only 1 count
		new_counts = []
		new_kmers = []
		for i in range(len(chunk_counts)):
			if chunk_counts[i] != 1:
				new_counts.append(chunk_counts[i])
				new_kmers.append(chunk_kmers[i])
		# Calculate size of remaining list and add this to list of distinct kmers for each chunk
		chunk_kmer_count = len(new_counts)
		kmer_count.append(chunk_kmer_count)
	# Calculate sum of different kmer counts from each chunk to get the total number of kmers for each run and write run accession number and kmer count to file
	kmer_count = sum(kmer_count)
	count_file.write(run+"\t"+kmer_count)
	subprocess.call(['rm', "/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/dump.fa"])
count_file.close()

