# Downloads fastq files from a list of SRA run accession numbers
# USAGE: python get-fastq.py [txt file containing list a SRA run accession numbers to download (one per line)]
import subprocess
import sys

with open(sys.argv[1]) as file:
	runs = file.readlines()
	runs = [line.strip('\n') for line in open(sys.argv[1])]

sixruns = list()

for run in runs:
	sixruns.append(run[0:6])

study = str(sys.argv[1][-13:-4])
newdir = "/media/nickeener/External_Drive/"+study
try:
	subprocess.call(['mkdir', newdir])
except:
	pass

for run,sixrun in zip(runs,sixruns):
	subprocess.call(['./get-fastq.sh', str(run), str(sixrun), newdir, study])