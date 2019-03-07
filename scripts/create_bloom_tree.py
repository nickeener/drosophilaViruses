# Automatically creates SBTs for reads from given SRA study and queries them to a the given reference
# USAGE: python create_bloom_tree.py [Study Accession Number] [0 for main drive || 1 for external drive]
import sys
import subprocess
import pandas as pd
import time

# Start timer
start_time = time.time()

# Set kmer size and drive to use
kmer = str(20)
drive = int(sys.argv[2])

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
print('ntcard --kmer='+kmer+' --threads=8 --pref='+study+' *.fastq.gz')
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
bf_size = str(bf_size)

# For each run, call makebf to create bloom filter from each run using the bit size calculated above
for run in runs:
	print('~/HowDeSBT/howdesbt makebf --k='+kmer+' --min=2 --bits='+bf_size+'K --threads=8 '+run+'_combined.fastq --out='+run+'.bf')
	subprocess.call(['./makebf.sh', study, run, kmer, bf_size])

# Run other howdesbt commands and query the index with the specified query
subprocess.call(['./cluster_build_query.sh', study, str(int(bf_size)*0.1)])

# Calculate and print total run time
time = time.time()-start_time
if time > 60 and time < 3600:
	hrs = 0
	mins = time//60
	secs = time/60
elif time >= 3600:
	hrs = time//3600
	mins = (time//3600)//60
	secs = (time//3600)%60
else:
	hrs = 0
	mins = 0
	secs = time

print('Total Time: '+str(hrs)+' Hours, '+str(mins)+' Minutes, '+str(secs)+' Seconds')