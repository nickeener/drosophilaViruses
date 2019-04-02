# Calculates number of Twyford kmers in each sequencing run that are present in at least 25% of other runs in a read set
# USAGE: python common_kmers.py [Study Accession Number] [Drive Type] [Query Name]
import sys
import subprocess
import pandas as pd
import pdb

# Set study accession number, kmer value and drive type to variables
study = sys.argv[1]
kmer = str(20)
drive = int(sys.argv[2])

# Set absolute path to query file to a variable
query = '/home/nickeener/projects/drosophilaViruses/viralSeqs/'+sys.argv[3]+'.fasta'

# Call ls_bf.sh to create file with name of each bloom filter on each line
subprocess.call(['./ls_bf.sh', study])

# Open file created by ls_bf.sh and read it line by line and store each line as element of a list
if drive == 0:
	with open('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/bf.txt') as file:
		bfs = file.read().splitlines()
	subprocess.call(['rm', '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/bf.txt'])
else:
	with open('/media/nickeener/External_Drive/'+study+'/bloomTrees/bf.txt') as file:
		bfs = file.read().splitlines()
	subprocess.call(['rm', '/media/nickeener/External_Drive/'+study+'/bloomTrees/bf.txt'])

# Create new list from bfs with the .bf removed from each to give the run accession number of each bloom filter
runs = []
for bf in bfs:
	runs.append(bf.strip('.bf'))

# Call querybf on each bloom filter in bfs
for bf, run in zip(bfs, runs):
	subprocess.call(['./querybf.sh', study, bf, query, run])

# For each run, read kmer file and store each kmer as an element of a list, then remove duplicates and store kmer list in a dictionary with its run number as its key
kmer_list = {}
for run in runs:
	if drive == 0:
		data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/'+run+'.kmer', header=None, sep=' ')
	else:
		data = pd.read_csv('/media/nickeener/External_Drive/'+study+'/bloomTrees/'+run+'.kmer', header=None, sep=' ')
	kmers = data[5].tolist()
	kmers = pd.unique(kmers).tolist()
	kmer_list[run] = kmers

# For each run, remove kmers that also appear in the Drosophila transcriptome (SRP111111), E. muscae transcriptome (SRP222222), and the wolbachia genome and create new dictionary to store unique kmers
data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP111111/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
background_kmers = data[5].tolist()
data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP222222/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
background_kmers.extend(data[5].tolist())
data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP333333/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
background_kmers.extend(data[5].tolist())
background_kmers = pd.unique(background_kmers).tolist()

unique_kmers = {}
for run in runs:
	unique_kmers[run] = list(set(kmer_list[run]) - set(background_kmers))

# Delete each individual kmer file
for run in runs:
	if drive == 0:
		kmer = '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/'+run+'.kmer'
	else:
		kmer = '/media/nickeener/External_Drive/'+study+'/bloomTrees/'+run+'.kmer'
	subprocess.call(['rm', kmer])

# Get list of all unique kmers from all runs
combined_kmers = []
for run in runs:
	kmers = unique_kmers[run]
	combined_kmers.extend(kmers)
combined_kmers = pd.unique(combined_kmers).tolist()

# Create dictionary where each kmer is a key and the corresponding value is the fraction of runs that the kmer is found in
kmer_freq = {}
for kmer in combined_kmers:
	kmer_freq[kmer] = 0
	for run in runs:
		if kmer in pd.unique(unique_kmers[run]).tolist():
			kmer_freq[kmer] += 1
for key in kmer_freq.keys():
	kmer_freq[key] = kmer_freq[key]/float(len(runs))

# Create new kmer dictionary using only kmers that are present in more than a quarter of the runs
new_kmer_freq = {}
for kmer, freq in kmer_freq.items():
	if freq >= 0.25:
		new_kmer_freq[kmer] = freq

# For each run, calculate the proportion of kmers that appear in at least 25% of all runs over the total number of viral kmers (8814)
kmers = new_kmer_freq.keys()
freqs = []
for run in runs:
	freqs.append(len(set(unique_kmers[run]).intersection(kmers))/float(8817))

if study == 'SRP133370':
	# Create data frame, sort my match proportion and write to csv
	sources = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/sources.txt', sep='\t', header=None)[1].tolist()
	data = {'Run Accession Number': runs, 'Match Proportion': freqs, "Source": sources}
	df = pd.DataFrame(data=data)
	df = df.sort_values(by=['Match Proportion'], ascending=False)
	df = df[['Run Accession Number', 'Match Proportion', 'Source']]
	df.to_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/kmer_proportions.txt', sep='\t', index=False, columns=['Run Accession Number', 'Match Proportion', 'Source'])
else:
	# Create data frame, sort my match proportion and write to csv
	data = {'Run Accession Number': runs, 'Match Proportion': freqs}
	df = pd.DataFrame(data=data)
	df = df.sort_values(by=['Match Proportion'], ascending=False)
	df = df[['Run Accession Number', 'Match Proportion']]
	if drive == 0:
		df.to_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomsTrees/kmer_proportions.txt', sep='\t', index=False, columns=['Run Accession Number', 'Match Proportion'])
	else:
		df.to_csv('/media/nickeener/External_Drive/'+study+'/bloomTrees/kmer_proportions.txt', sep='\t', index=False, columns=['Run Accession Number', 'Match Proportion'])


'''# Get list of all distinct kmers that in Exposed Whole Fly or Infected Cadaver runs and another list of all distinct kmers in Unexposed Whole Fly runs
sources = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/sources.txt', sep='\t', header=None)
exposed_runs = sources.loc[(sources[1] == 'Exposed Whole Fly') | (sources[1] == "Infected Cadaver")]
exposed_runs = exposed_runs[0].tolist() 
unexposed_runs = sources.loc[sources[1] == 'Unexposed Whole Fly']
unexposed_runs = unexposed_runs[0].tolist()

exposed_kmers_set = []
unexposed_kmers_set = []
for run in exposed_runs:
	exposed_kmers_set.extend(set(unique_kmers[run]))
exposed_kmers_set = pd.unique(exposed_kmers_set).tolist()
for run in unexposed_runs:
	unexposed_kmers_set.extend(set(unique_kmers[run]))
unexposed_kmers_set = pd.unique(unexposed_kmers_set).tolist()

# Create dictionary from each list with the kmer as the key and the number of runs it occurs in as the value
exposed_kmer_freq = {}
for kmer in exposed_kmers_set:
	exposed_kmer_freq[kmer] = 0
	for run in exposed_runs:
		if kmer in unique_kmers[run]:
			exposed_kmer_freq[kmer] += 1
for key in exposed_kmer_freq.keys():
	exposed_kmer_freq[key] = exposed_kmer_freq[key]/float(len(exposed_runs))

unexposed_kmer_freq = {}
for kmer in unexposed_kmers_set:
	unexposed_kmer_freq[kmer] = 0
	for run in unexposed_runs:
		if kmer in unique_kmers[run]:
			unexposed_kmer_freq[kmer] += 1
for key in unexposed_kmer_freq.keys():
	unexposed_kmer_freq[key] = unexposed_kmer_freq[key]/float(len(unexposed_runs))

# Create new dictionaries using only kmers that are present in more than 25% of runs
new_exposed_kmer_freq = {}
for kmer,freq in exposed_kmer_freq.items():
	if freq >= 0.25:
		new_exposed_kmer_freq[kmer] = freq

new_unexposed_kmer_freq = {}
for kmer,freq in unexposed_kmer_freq.items():
	if freq >= 0.25:
		new_unexposed_kmer_freq[kmer] = freq

# For each exposed/unexposed run, calculate the proportion of kmers that appear in at least 25% of all runs over the total number of viral kmers (8814)
kmers = new_exposed_kmer_freq.keys()
exposed_freqs = []
for run in exposed_runs:
	exposed_freqs.append(len(set(unique_kmers[run]).intersection(kmers))/float(8814))

#kmers = new_unexposed_kmer_freq.keys()
unexposed_freqs = []
for run in unexposed_runs:
	unexposed_freqs.append(len(set(unique_kmers[run]).intersection(kmers))/float(8814))

# Create data frame, sort my match proportion and write to csv (for both exposed and unexposed runs)
data = {'Run Accession Number': exposed_runs, 'Match Proportion': exposed_freqs}
exposed_df = pd.DataFrame(data=data)
exposed_df = exposed_df.sort_values(by=['Match Proportion'], ascending=False)
exposed_df = exposed_df[['Run Accession Number', 'Match Proportion']]

data = {'Run Accession Number': unexposed_runs, 'Match Proportion': unexposed_freqs}
unexposed_df = pd.DataFrame(data=data)
unexposed_df = unexposed_df.sort_values(by=['Match Proportion'], ascending=False)
unexposed_df = unexposed_df[['Run Accession Number', 'Match Proportion']]'''
