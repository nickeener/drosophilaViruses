# Creates file with the number of reads in each accession for a set of reads
# USAGE: python readcounts.py [Study Accession Number] [Drive Type]

# Import libraries
import sys
import subprocess
import pdb

# Set study accession to and path to mapping directory to variables
study = sys.argv[1]
if sys.argv[2] == 0:
	mapdir = '/home/nickeener/projects/drosophilaViruses/mapping/'+study
else:
	mapdir = '/media/nickeener/External_Drive'

# Get list of all read files (one per pair) and create list with the absolute path to each file
ls = subprocess.Popen(['ls', mapdir+'/*_1.fastq.gz'], stdout=subprocess.PIPE)
runs = ls .stdout.read()

for run in runs:
	subprocess.call(['./readcounts.sh', study, run])