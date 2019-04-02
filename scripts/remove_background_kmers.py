# Uses querybf on simmed Drosophila reads to find kmers shared between Drosophila and the query then uses querybf again to find twyford kmers in actual sequencing data
# Background kmers are then removed from data and the number of matches recalculated and the output of howdesbt query modified to reflect the new kmer match count
# USAGE: python remove_background_kmers.py [Study Accession Number] [Drive Type] [Query Name]

import sys
import subprocess
import pandas as pd
import pdb

# Set study accession number and kmer value to a variable
study = sys.argv[1]
kmer = str(20)

# Set absolute path to query file to a variable
query = '/home/nickeener/projects/drosophilaViruses/viralSeqs/'+sys.argv[3]+'.fasta'

# Set drive type (0 for main drive, 1 for external drive)
drive = int(sys.argv[2])

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
	subprocess.call(['rm', '/media/nickeener/External_Drive/'+study+'bloomTrees/bf.txt'])

# Create new list from bfs with the .bf removed from each to give the run accession number of each bloom filter
runs = []
for bf in bfs:
	runs.append(bf.strip('.bf'))

# Call querybf on each bloom filter in bfs
for bf, run in zip(bfs, runs):
	subprocess.call(['./querybf.sh', study, bf, query, run]) 

# Only run the following line when creating new background kmer lists
# Combine kmer files to create file containing all kmers
'''subprocess.call(['./combine_bfs.sh', study, kmer])'''

# For creating new background kmer lists, all subsequent lines should be commented out
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

# For each run, remove kmers that also appear in Drosophila (SRP111111) and E, muscae (SRP222222) transcriptomes and the Wolbachia (SRP333333), Sigma Virus (SRP444444),
# Drosophila C Virus (SRP555555), and Nora Virus (SRP666666) genomes and then and create new dictionary to store unique kmers
data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP111111/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
background_kmers = data[5].tolist()
data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP222222/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
background_kmers.extend(data[5].tolist())
data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP333333/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
background_kmers.extend(data[5].tolist())
try:
	data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP444444/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
	background_kmers.extend(data[5].tolist())
except:
	pass
try:
	data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP555555/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
	background_kmers.extend(data[5].tolist())
except:
	pass
try:
	data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/SRP666666/bloomTrees/combined_k'+kmer+'.kmer', header=None, sep=' ')
	background_kmers.extend(data[5].tolist())
except:
	pass
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

# Modify the .dat file to reflect the new kmer count
# Read .dat file
if drive == 0:
	data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/queries_'+sys.argv[3]+'.dat', header=None, sep=' ')
else:
	data = pd.read_csv('/media/nickeener/External_Drive/'+study+'/bloomTrees/queries_'+sys.argv[3]+'.dat', header=None, sep=' ')

# Creates lists with each match fraction and its corresponding run accession number
match_fraction = data[1].tolist()
match_fraction = match_fraction[1:len(match_fraction)] # Removes first entry which is from header line
match_runs = data[0].tolist()
match_runs = match_runs[1:len(match_runs)] # Same deal as above

# Extract out denominator from first fraction (should be the same for all)
total_kmers = ''
for i in range(len(match_fraction[0])):
	if match_fraction[0][i] != '/':
		pass
	else:
		j = i+1
		break
for i in match_fraction[0][j:]:
	total_kmers += i

# Create new match fraction list using the number of unique kmers and calculate the decimal percentage for each and store in a list
new_match_fractions = []
new_match_decimals = []
for match_run in match_runs:
	new_match_fractions.append(str(len(unique_kmers[match_run]))+'/'+total_kmers)
	decimal = len(unique_kmers[match_run])/float(total_kmers)
	new_match_decimals.append(str(decimal))

# Get source data (if study is SRP133370)
if study == 'SRP133370':
	sources = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/sources.txt', sep='\t')['Source'].tolist()

# Create dataframe from match_runs, new_match_fractions, and new_match_decimals and write to csv
data = {'Run Accession Number': match_runs, 'Match Fraction': new_match_fractions, 'Decimal Fraction': new_match_decimals}
df = pd.DataFrame(data=data)
df = df.sort_values(by=['Run Accession Number'])
match_runs = df['Run Accession Number']
new_match_fractions = df['Match Fraction']
new_match_decimals = df['Decimal Fraction']
if study == 'SRP133370':
	data = {'Run Accession Number': match_runs, 'Match Fraction': new_match_fractions, 'Decimal Fraction': new_match_decimals, 'Source': sources}
	df = pd.DataFrame(data=data) 
	df = df[['Run Accession Number', 'Match Fraction', 'Decimal Fraction', 'Source']]
else:
	data = {'Run Accession Number': match_runs, 'Match Fraction': new_match_fractions, 'Decimal Fraction': new_match_decimals}
	df = pd.DataFrame(data=data) 
	df = df[['Run Accession Number', 'Match Fraction', 'Decimal Fraction']]
df = df.sort_values(by=['Decimal Fraction'], ascending=False)
if drive == 0 and study == 'SRP133370':
	df.to_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/corrected_query_'+sys.argv[3]+'.dat', sep='\t', index=False, columns=['Run Accession Number', 'Match Fraction', 'Decimal Fraction', 'Source'])
elif drive == 0 and study != 'SRP133370':
	df.to_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees/corrected_query_'+sys.argv[3]+'.dat', sep='\t', index=False, columns=['Run Accession Number', 'Match Fraction', 'Decimal Fraction'])
else: 
	df.to_csv('/media/nickeener/External_Drive/'+study+'/bloomTrees/corrected_query_'+sys.argv[3]+'.dat', sep='\t', index=False, columns=['Run Accession Number', 'Match Fraction', 'Decimal Fraction'])

