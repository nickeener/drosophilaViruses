# Reads input files (consisting of the mate pair file names listed on separate lines), converts
# this information into a format that can be inputed into bowtie2 and then calls bowtie2
# to map the read files to the specified reference
# USAGE: python bowtie2.py <Study Accession Number> <Reference Index>
	# Study Accession Number is the folder name where the read files are located
	# Reference Index is the folder name where the index files for the reference genome 
	# are located

# Make sure read file list file is formatted correctly! (No spaces at end of each)

import sys
import subprocess

# Read file, store as list and remove newline characters
pair1file = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/pair1_trimmed.txt"
with open(pair1file) as file:
	pair1 = file.readlines()
	pair1 = [line.strip('\n') for line in open(pair1file)]

# Create list with the absolute path to each file
pairpath = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/trimmed/"
newpair1 = []
for i in pair1:
	newpair1.append(pairpath+i)

# Convert previous list into a comma separated string
pair1commalist= ""
for i in newpair1:
	pair1commalist += i+","
pair1commalist = pair1commalist[:-1]

# Read file, store as list and remove newline characters
pair2file = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/pair2_trimmed.txt"
with open(pair2file) as file:
	pair2 = file.readlines()
	pair2 = [line.strip('\n') for line in open(pair2file)]

# Create list with the absolute path to each file
newpair2 = []
for i in pair2:
	newpair2.append(pairpath+i)


# Convert previous list into a comma separated string
pair2commalist= ""
for i in newpair2:
	pair2commalist += i+","
pair2commalist = pair2commalist[:-1]

# Store absolute path to reference genome index
reference = sys.argv[2]
index = "/home/nickeener/projects/drosophilaViruses/mapping/indexes/"+reference+"/bowtie2/"+reference

# Create new folder to store results in
newdir = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/bowtie"
try:
	subprocess.call(['mkdir', newdir])
except:
	pass

# Call bowtie2
subprocess.call(['./bowtie2.sh', index, pair1commalist, pair2commalist, newdir, sys.argv[2]])

# Debugging
'''print(index)
print("\n")
print(pair1commalist)
print("\n")
print(pair2commalist)
print("\n")
print(newdir)'''

#print("bowtie2 -t --very-sensitive-local -p 8 -x "+index+" -1 "+pair1commalist+" -2 "+pair2commalist+" -S "+newdir+"/"+sys.argv[2]+".sam")

