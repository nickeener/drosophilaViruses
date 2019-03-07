# Usage: python readcounts.py [Study Acession Number] [Drive Type (0 for main drive or 1 for external drive)]
import sys
import subprocess
import pdb

# Assign study accession number to variable and set drive type
study = sys.argv[1]
drive = int(sys.argv[2])

# Get list of all files/directories in the mapping directory of interest
if drive == 0:
	ls = subprocess.Popen(['ls', '/home/nickeener/projects/drosophilaViruses/mapping/'+study], stdout=subprocess.PIPE)
else:
	ls = subprocess.Popen(['ls', '/media/nickeener/External_Drive/'+study], stdout=subprocess.PIPE)
output = ls.stdout.read()

# Convert output string into list containing only the read files
reads = []
for i in range(len(output)):
	if output[i] == 'S' or output[i] == 'E' and output[i+1] == 'R':
		reads.append(output[i:i+21])

# Elimnate every other element from previous list (each run has a double because of forward/reverse reads)
runs = []
for i in range(len(reads)):
	if i%2 == 1:
		runs.append(reads[i])

# Call readcounts.sh on each run to print out the number of reads in each run to stdout
for run in runs:
	print(run)
	output = subprocess.Popen(['./readcounts.sh', study, run], stdout=subprocess.PIPE)
	linecount = output.stdout.read()
	print(int(linecount)/4)