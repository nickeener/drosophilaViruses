# Automatically creates SBTs for reads from given SRA study and queries them to a the given reference
# USAGE: python create_bloom_tree.py [Study Accession Number] [0 for main drive || 1 for external drive]
import sys
import subprocess
#import pandas as pd
import time
import datetime
import numpy
import math
import pdb

# Start timer
start_time = time.time()

# Set kmer size and drive to use
kmer = str(20)
drive = int(sys.argv[2])
bf_size = str(878967)

# Get list of all files/directories in the mapping directory of interest
study = sys.argv[1]
if drive == 0:
	ls = subprocess.Popen(['ls', '/home/nickeener/projects/drosophilaViruses/mapping/'+study], stdout=subprocess.PIPE)
else:
	ls = subprocess.Popen(['ls', '/media/nickeener/External_Drive/'+study], stdout=subprocess.PIPE)
output = ls.stdout.read()

# Convert output string into list containing only the SRA run accession number of each file
reads = []
for i in range(len(output)):
	if output[i] == 'S' or output[i] == 'E' and output[i+1] == 'R':
		reads.append(output[i:i+10])

# Elimnate every other element from previous list (each run has a double because of forward/reverse reads)
runs = []
for i in range(len(reads)):
	if i%2 == 1:
		runs.append(reads[i])

# Create bloomTree directory in correct mapping directory
if drive == 0:
	newdir = '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees'
else:
	newdir = '/media/nickeener/External_Drive/'+study+'/bloomTrees'
subprocess.call(['mkdir', newdir])

# Calculate kmer frequencies using ntcard and use to calculate appropriate bloom filter size
'''print('ntcard --kmer='+kmer+' --threads=8 --pref='+study+' *.fastq.gz')
subprocess.call(['./ntcard.sh', kmer, study])
if drive == 0:
	data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'+study+'_k'+kmer+'.hist', header=None, sep='\t', names=['1', '2'])
else:
	data = pd.read_csv('/media/nickeener/External_Drive/'+study+'/'+study+'_k'+kmer+'.hist', header=None, sep='\t', names=['1', '2'])
bf_size = data['2'][1]-data['2'][2]
if drive == 0:
	subprocess.call(['rm', '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'+study+'_k'+kmer+'.hist'])
else:
	subprocess.call(['rm', '/media/nickeener/External_Drive/'+study+'/'+study+'_k'+kmer+'.hist'])
bf_size = (bf_size+int(bf_size*.05))/1000 # Add a small portion to make a slight overstimation and divide by 1000 to get it in K format
bf_size = str(bf_size)'''


# For each run, call makebf to create bloom filter from each run using the bit size calculated above
'''for run in runs:
	print('~/tools/HowDeSBT/howdesbt makebf --k='+kmer+' --min=2 --bits='+bf_size+'K --threads=8 '+run+'_combined.fastq --out='+run+'.bf')
	subprocess.call(['./makebf.sh', study, run, kmer, bf_size])'''

# Run other howdesbt commands and query the index with the specified query
#subprocess.call(['./cluster_build_query.sh', study, str(int(bf_size)*0.1)])

# Calculate and print total run time
time = time.time()-start_time
newtime = str(datetime.timedelta(seconds=time))
#print('Total Runtime: '+newtime)

# Read output file and use regression algorithm to find average of "empty" (no virus) runs

# Open file and remove trailing newlines
with open(newdir+'/queries_drosophilaViruses.dat') as file:
	lines = file.readlines()
for i in range(len(lines)):
	lines[i] = lines[i].rstrip()

# Create a list where each element is a list with two elements, a string with the viral sequence name and a dictionary 
# with run accession numbers as keys and the corresponding kmer match proportions as their values
matches = []
i = 0
count = 0
while i < len(lines):
	if lines[i][0] == '*':
		matches.append([lines[i], {}])
		j = i
		j += 1
		while j < len(lines) and lines[j][0] != '*':
			j += 1
	for k in range(i+1, j):
		matches[count][1][lines[k].split()[0]] = lines[k].split()[2]
	i = j
	count += 1

# Use regression algorithm to calculate average of uninfected runs
averages = []
for match in matches:
	values = []
	values.append(match[1].values())
	log_values = []
	for value in values[0]:
		log_values.append(math.log10(float(value)))
	averages.append(numpy.mean(log_values))

average = numpy.mean(averages)
std = numpy.std(averages)



new_average = average + 0.001
new_averages = []
while round(average, 5) != round(new_average, 5):
	average = new_average
	for avg in averages:
		if avg <= average+(2*std) and avg >= math.log10(0.05):
			new_averages.append(avg)
	averages = new_averages
	new_average = numpy.mean(averages)
	std = numpy.std(averages)
	new_averages = []

for average in averages:
	print(average)
