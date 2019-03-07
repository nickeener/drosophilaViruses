import sys
import subprocess

# Read file, store as list and remove newline characters
pair1file = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/pair1_trimmed50_noqual.txt"
with open(pair1file) as file:
	pair1 = file.readlines()
	pair1 = [line.strip('\n') for line in open(pair1file)]

# Create list with the absolute path to each file
pairpath = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/trimmed50_noqual/"
newpair1 = []
for i in pair1:
	newpair1.append(pairpath+i)

# Read file, store as list and remove newline characters
pair2file = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/pair2_trimmed50_noqual.txt"
with open(pair2file) as file:
	pair2 = file.readlines()
	pair2 = [line.strip('\n') for line in open(pair2file)]

# Create list with the absolute path to each file
newpair2 = []
for i in pair2:
	newpair2.append(pairpath+i)

run_acc = []
for pair in pair1:
	run_acc.append(pair[:10])


for pair1, pair2, run in zip(newpair1, newpair2, run_acc):
	subprocess.call(['./call_bbmerge.sh', pair1, pair2, run])
