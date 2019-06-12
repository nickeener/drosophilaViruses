# Downloads fastq files from a list of SRA run accession numbers
# USAGE: python get-fastq.py [txt file containing list a SRA run accession numbers to download (one per line)] [Drive Type]
import subprocess
import sys
import pdb

with open(sys.argv[1]) as file:
	runs = file.readlines()
	runs = [line.strip('\n') for line in open(sys.argv[1])]

sixruns = list()

for run in runs:
	sixruns.append(run[0:6])

study = str(sys.argv[1][-13:-4])

if int(sys.argv[2]) == 0:
	newdir= "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1][:-4]+"/"
else:
	newdir = "/media/nickeener/External_Drive/"+sys.argv[1][:-4]+"/"
try:
	subprocess.call(['mkdir', newdir])
except:
	pass

for run,sixrun in zip(runs,sixruns):
	print("Fetching...."+run)
	subprocess.call(['./get-fastq.sh', str(run), str(sixrun), newdir, study])