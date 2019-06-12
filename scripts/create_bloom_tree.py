# Automatically creates SBTs for reads from given SRA study and queries them to a the given reference
# USAGE: python create_bloom_tree.py [Study Accession Number] [0 for main drive || 1 for external drive] [Query Name]
import sys
import subprocess
#import pandas as pd
import time
import datetime
import numpy
import math
import pdb

# Start timer
start_time = time.time()

# Set kmer size and drive to use
kmer = str(20)
drive = int(sys.argv[2])
bf_size = str(878967)

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
'''print('ntcard --kmer='+kmer+' --threads=8 --pref='+study+' *.fastq.gz')
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
bf_size = str(bf_size)'''


# For each run, call makebf to create bloom filter from each run using the bit size calculated above
'''for run in runs:
	print('~/tools/HowDeSBT/howdesbt makebf --k='+kmer+' --min=2 --bits='+bf_size+'K --threads=8 '+run+'_combined.fastq --out='+run+'.bf')
	subprocess.call(['./makebf.sh', study, run, kmer, bf_size])'''

# Run other howdesbt commands and query the index with the specified query
#subprocess.call(['./SBT_cluster_build.sh', study, str(int(bf_size)*0.1)])
#subprocess.call(['SBT_query', study, 'drosophilaViruses'])

# Calculate and print total run time
time = time.time()-start_time
newtime = str(datetime.timedelta(seconds=time))
print('Total Runtime: '+newtime)

# For each virus, query the SBT with a set of 100 random permutations of the sequence and calculate a baseline value
# for each virus
# AND create a list where each element is a list with two elements, a string with the viral sequence name and a dictionary 
# with run accession numbers as keys and the corresponding kmer match proportions as their values
seqdir = '/home/nickeener/projects/drosophilaViruses/viralSeqs/'
ls = subprocess.Popen(['ls', seqdir+'permutations'], stdout=subprocess.PIPE)
perms = ls.stdout.read()
perms = perms.splitlines()

'''for perm in perms:
	subprocess.call(['./SBT_perm_query.sh', study, 'permutations/'+perm, perm])'''

cutoffs = []
matches = []

for perm,j in zip(perms,range(len(perms))):
	marker = 0
	with open(newdir+'/queries_'+perm+'.dat') as file:
		lines = file.readlines()
	for i in range(len(lines)):
		lines[i] = lines[i].rstrip()
	values = []
	for line in lines:
		if list(line)[0] == '*':
			pass
		else:
			values.append(math.log10(float(line.split()[2])))
	'''mean = numpy.mean(values)
	std = numpy.std(values)'''

	'''# Regression algorithm
	new_mean = mean + 0.00000001
	new_values = []
	while round(mean, 10) != round(new_mean, 10):
		mean = new_mean
		for value in values:
			if value <= mean+(3*std):
				new_values.append(value)
		values = new_values
		new_mean = numpy.mean(values)
		std = numpy.std(values)
		new_values = []
	cutoffs.append([lines[0], 10**(mean+(3*std))])'''
	cutoffs.append([lines[0], 10**(numpy.quantile(values, 0.99))])

	# For each virus, write to a file all the values that passed this step
	with open(newdir+'/distributions/'+lines[0].split()[0][1:]+'.txt', 'w') as file:
		for value in values:
			file.write(str(value)+'\n')


# Create a list where each element is a list with two elements, a string with the viral sequence name and a dictionary 
# with run accession numbers as keys and the corresponding kmer match proportions as their values

# Open file and remove trailing newlines
with open(newdir+'/queries_twyford_full.dat') as file:
	lines = file.readlines()
for i in range(len(lines)):
	lines[i] = lines[i].rstrip()

matches = []
i = 0
count = 0
while i < len(lines):
	if lines[i][0] == '*':
		matches.append([lines[i], {}])
		j = i
		j += 1
		while j < len(lines) and lines[j][0] != '*':
			j += 1
	for k in range(i+1, j):
		matches[count][1][lines[k].split()[0]] = lines[k].split()[2]
	i = j
	count += 1

# Create list (similar format to matches) where each element is a list with element 1 as the virus name and element 2 
# as a dictonary with the SRA run accession numbers as keys and their kmer match proportion for the particular virus
# as the values (if that value exceeds the cutoff value)
hits = []
for match,i in zip(matches, range(len(matches))):
	if int(match[0][-2:]) >= 10:
		hits.append([match[0][:-2], {}])
	else:
		hits.append([match[0][:-1], {}])
	for key in match[1].keys():
		for cutoff in cutoffs:
			if float(match[1][key]) >= cutoff[1]:
				hits[i][1][key] = match[1][key]

# Create a dictonary with only the top kmer match proportion for each accession
top_hits = {}
for hit in hits:
	for key in hit[1].keys():
		if key in top_hits.keys():
			if hit[1][key] > top_hits[key]:
				top_hits[key] = hit[1][key]
			else:
				pass
		else:
			top_hits[key] = hit[1][key]

# Create matrix of significant hits
sig_hits = []
for hit,i  in zip(hits, range(len(hits))):
	sig_hits.append([hit[0], {}])
	for key in hit[1].keys():
		if hit[1][key] == top_hits[key]:
			sig_hits[i][1][key] = hit[1][key]

# Write significant hits to file (in decreasing order)
subprocess.call(['rm', newdir+'/new_queries_twyford_full.dat'])
with open(newdir+'/new_queries_twyford_full.dat', 'w') as file:
	for hit in sig_hits:
		file.write(hit[0]+''+str(len(hit[1]))+'\n')
		sort = hit[1].values()
		new_sort = []
		for value in sort:
			new_sort.append(float(value))
		new_sort.sort(reverse=True)
		for value in new_sort:
			for key in hit[1].keys():
				if float(hit[1][key]) == value:
					file.write(key+' '+str(value)+'\n')
		file.write('\n')
	file.write('\nEND SIGNIFICANT HITS\n\n')
	for line in lines:
		file.write(line+'\n')









'''# Read output file and use regression algorithm to find average of "empty" (no virus) runs (OLD WAY)

# Open file and remove trailing newlines
with open(newdir+'/queries_drosophilaViruses.dat') as file:
	lines = file.readlines()
for i in range(len(lines)):
	lines[i] = lines[i].rstrip()

# Create a list where each element is a list with two elements, a string with the viral sequence name and a dictionary 
# with run accession numbers as keys and the corresponding kmer match proportions as their values
matches = []
i = 0
count = 0
while i < len(lines):
	if lines[i][0] == '*':
		matches.append([lines[i], {}])
		j = i
		j += 1
		while j < len(lines) and lines[j][0] != '*':
			j += 1
	for k in range(i+1, j):
		matches[count][1][lines[k].split()[0]] = lines[k].split()[2]
	i = j
	count += 1

# Use regression algorithm to calculate average of uninfected runs

# Extract all values and store their log10 values in a list, then calculate initial mean and standard deviation
values = []
for match in matches:
	for key in match[1]:
		values.append(math.log10((float(match[1][key]))))
average = numpy.mean(values)
std = numpy.std(values)

# Calculate average and standard deviation of all values that are not likely to contain viral sequeneces by excluding 
# values greater than 3 standard deviations away from the mean and then recalculating the average and standard deviation 
# and repeating until the average stops changing.
new_average = average+.00001
new_values = []
final_values = []
while round(average, 10) != round(new_average, 10):
	average = new_average
	for value in values:
		if value <= average+(3*std) and value <= math.log10(0.5):
			new_values.append(value)
	values = new_values
	final_values = new_values
	new_average = numpy.mean(values)
	std = numpy.std(values)
	new_values = []

cutoff = 10**(average+(3*std))

# Update the .dat file so it has a list of the runs that had significant matches at the top, reporting only the top hit for
# each run
subprocess.call(['rm', newdir+'/new_queries_drosophilaViruses.dat'])
with open(newdir+'/new_queries_drosophilaViruses.dat', 'w') as file:
	file.write('SIGNIFICANT HITS\nUsing a cutoff of: '+str(round(cutoff, 6))+'\n\n')
	all_hits = []
	# Create list of all hits that have a kmer match proportion higher than the cutoff value
	for match in matches:
		if int(match[0][-2:]) >= 10:
			hits = [[match[0][:-2]], {}]
		else:
			hits = [[match[0][:-1]], {}]
		for value in match[1].values():
			if float(value) >= cutoff:
				for key in match[1].keys():
					if match[1][key] == value:
						hits[1][key] = value
		all_hits.append(hits)
	# If a run has a significant hit for more than one virus, key only the highest scoring hit
	sig_hits = []
	all_keys = []
	for hit in all_hits:
		if len(hit[1]) == 0:
			pass
		else:
			keys = hit[1].keys()
		try:
			for key in keys:
				all_keys.append(key)
		except:
			pass
	uniq_keys = set(all_keys) # Gets all unique run accession numbers
	# Finds max value of each run accession number and adds only that accession : value pair to the sig_hits list
	for hit in all_hits:
		sig_hits.append([hit[0], {}])
	for key in uniq_keys:
		max_value = 0
		for hit in all_hits:
			try:
				if hit[1][key] > max_value:
					max_value = hit[1][key]
			except:
				pass
		for i in range(len(all_hits)):
			for key1 in all_hits[i][1].keys():
				try:
					if all_hits[i][1][key] == max_value:
						sig_hits[i][1][key] = max_value
				except:
					pass
	# Write significant hits to file (in decreasing order)
	for hit in sig_hits:
		file.write(hit[0][0]+' '+str(len(hit[1]))+'\n')
		sort = hit[1].values()
		new_sort = []
		for value in sort:
			new_sort.append(float(value))
		new_sort.sort(reverse=True)
		for value in new_sort:
			for key in hit[1].keys():
				if float(hit[1][key]) == value:
					file.write(key+' '+str(value)+'\n')
		file.write('\n')
	file.write('\nEND SIGNIFICANT HITS\n\n')
	for line in lines:
		file.write(line+'\n')'''