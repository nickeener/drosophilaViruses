# Reads input files (consisting of the mate pair file names listed on separate lines), converts
# this information into a format that can be inputed into bowtie2 and then calls bowtie2
# to map the read files to the specified reference
# USAGE: python bowtie2.py <Study Accession Number> <Reference Index> <Drive Type>
	# Study Accession Number is the folder name where the read files are located
	# Reference Index is the folder name where the index files for the reference genome 
	# are located

# Make sure read file list file is formatted correctly! (No spaces at end of each)

# Import libraries
import sys
import subprocess
import pdb

# Set study accession number and drive type to variable
study = sys.argv[1]
drive = int(sys.argv[3])
print("0")
# Read file, store as list and remove newline characters
if drive == 0:
	pair1file = "/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/pair1.txt"
else:
	pair1file = "/media/nickeener/External_Drive/"+study+"/pair1.txt"
with open(pair1file) as file:
	pair1 = file.readlines()
	pair1 = [line.strip('\n') for line in open(pair1file)]

# Create list with the absolute path to each file and sort list
if drive == 0:
	pairpath = "/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/"
else:
	pairpath = "/media/nickeener/External_Drive/"+study+"/"
newpair1 = []
for i in pair1:
	newpair1.append(pairpath+i)
newpair1.sort()

# Convert previous list into a comma separated string
pair1commalist= ""
for i in newpair1:
	pair1commalist += i+","
pair1commalist = pair1commalist[:-1]

# Read file, store as list and remove newline characters
if drive == 0:
	pair2file = "/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/pair2.txt"
else:
	pair2file = "/media/nickeener/External_Drive/"+study+"/pair2.txt"
with open(pair2file) as file:
	pair2 = file.readlines()
	pair2 = [line.strip('\n') for line in open(pair2file)]

# Create list with the absolute path to each file and sort list
newpair2 = []
for i in pair2:
	newpair2.append(pairpath+i)
newpair2.sort()

# Convert previous list into a comma separated string
pair2commalist= ""
for i in newpair2:
	pair2commalist += i+","
pair2commalist = pair2commalist[:-1]

# Store absolute path to reference genome index
reference = sys.argv[2]
index = "/home/nickeener/projects/drosophilaViruses/mapping/indexes/bowtie2/"+reference

# Create new folder to store results in
if drive == 0:
	newdir = "/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/bowtie"
else:
	newdir = "/media/nickeener/External_Drive/"+study+"/bowtie"
try:
	subprocess.call(['mkdir', newdir])
except:
	pass

# Call bowtie2 on all runs in a study at once
subprocess.call(['./bowtie2.sh', index, pair1commalist, pair2commalist, newdir, sys.argv[2]])
	
'''# Run the following lines only if aligning each fastq file individually	
if drive == 0:
	newdir = "/home/nickeener/projects/drosophilaViruses/mapping/"+study+"/bowtie/individual"
else:
	newdir = "/media/nickeener/External_Drive/"+study+"/bowtie/individual"
try:
	subprocess.call(['mkdir', newdir])
except:
	pass


# Call bowtie2 on each run individually
for pair1,pair2 in zip(newpair1,newpair2):
	if drive == 0:
		print("Mapping "+study+" - "+pair1[61:71])
		subprocess.call(['./bowtie2.sh', index, pair1, pair2, newdir, pair1[61:71], reference])
	else:
		print("Mapping "+pair1[42:52])
		subprocess.call(['./bowtie2.sh', index, pair1, pair2, newdir, pair1[42:52], reference])'''


# Debugging
'''print(index)
print("\n")
print(pair1commalist)
print("\n")
print(pair2commalist)
print("\n")
print(newdir)'''

#print("bowtie2 -t --very-sensitive-local -p 8 -x "+index+" -1 "+pair1commalist+" -2 "+pair2commalist+" -S "+newdir+"/"+sys.argv[2]+".sam")

