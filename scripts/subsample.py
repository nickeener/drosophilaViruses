# USAGE: python subsample.py [Study Accession Number] [Drive Type]
# Creates subsample of all reads in a study
import sys
import subprocess

# Set study accession number, drive type and samplerate to variables
study = sys.argv[1]
drive = int(sys.argv[2])
samplerate = [str(0.05), str(0.15), str(0.2), str(0.25), str(0.4)]

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
		reads.append(output[i:i+21])

# Plave pair1 variables in one list and pair2 variables in another
pair1 = []
pair2 = []
for i in range(len(reads)):
	if i%2 == 0:
		pair1.append(reads[i])
	else:
		pair2.append(reads[i])

for rate in samplerate:
# Create new directory to put subsampled reads in for each samplerate
	if drive == 0:
		newdir = '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/subsample'+rate
	else:
		newdir = '/media/nickeener/External_Drive/'+study+'/subsample'+rate
	subprocess.call(['mkdir', newdir])


# Call subsample.sh which uses reformat.sh (bbmap) to create subsamples of reads in a study for each samplerate
for rate in samplerate:
	for one,two in zip(pair1,pair2):
		subprocess.call(['./subsample.sh', study, one, two, rate])
		print("Done with "+one[0:10])