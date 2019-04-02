# Usage python large_create_bloom_tree.py [Study Accession Number] [Run Accession Number] [Drive Type]

import sys
import subprocess
#import pandas as pd
import shlex
import pdb

# Assign study and run accession numbers, kmer value and drive type to variables
study = sys.argv[1]
run = sys.argv[2]
kmer = str(20)
drive = int(sys.argv[3])
#bf_size = str(834090)

# Create bloomTree directory in correct mapping directory and store path to read directory in variable
if drive == 0:
	newdir = '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees'
	rundir = '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'
else:
	newdir = '/media/nickeener/External_Drive/'+study+'/bloomTrees'
	rundir = '/media/nickeener/External_Drive/'+study+'/'
subprocess.call(['mkdir', newdir])

# Calculate kmer frequencies of the large read file using ntcard and use to calculate appropriate bloom filter size
command = shlex.split('ntcard --kmer='+kmer+' --threads=8 --pref='+run+' '+rundir+run+'_*.fastq.gz')
print(command)
subprocess.run(command, shell=True)
#subprocess.call(shlex.split('ntcard --kmer='+kmer+' --threads=8 --pref='+run+ ' '+rundir+run+'_*.fastq.gz'))
#pdb.set_trace()