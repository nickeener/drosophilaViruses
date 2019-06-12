# Python script that extracts run accession number and read count data from bowtie2 log files
# Creates tabular file with run accession numbers in one column, the number of mapped reads in another column, and 
# and the total number of reads in another column
# USAGE python extract_read_counts.py [Study Accession Number] [Drive Type] > [File Name]

# Libraries
import sys
import subprocess
import pandas as pd
import pdb

# Set drive type
drive = int(sys.argv[2])

# Assign path of directory containing .log files to a variable, then get names of all .log files and store in list
path = "/home/nickeener/projects/drosophilaViruses/mapping/"+sys.argv[1]+"/bowtie/individual/"
out = subprocess.Popen(['ls', path], stdout=subprocess.PIPE)
logs = out.stdout.readlines()
newlogs = []
for log in logs:
	newlogs.append(str(log.rstrip())) # Removing trailing \n's
if 'readcounts.txt' in newlogs:
	del newlogs[newlogs.index('readcounts.txt')] # Removes old readcount file if it exists
logs = []
for i in range(0, len(newlogs), 2): 
	logs.append(str(newlogs[i])) # Removes SAM files from list

# Read each .log file in the specified directory and store the run accession number, mapped read count, and total reads count
# in separate lists
run_acc = []
mapped_count = []
total_count = []
for log in logs:
	run_acc.append(log[:-4])
	with open(path+log) as file:
		output = file.readlines()
		mapped_count.append(output[7].split(' ')[4])
		total_count.append(output[4].split(' ')[0])

# Create tabular file containing all the data
data = {'Run Accession Number': run_acc, 'Mapped Read Count': mapped_count, 'Total Read Count': total_count}
df = pd.DataFrame(data=data)
df.to_csv(path+'readcounts.txt', sep = '\t', index=False)