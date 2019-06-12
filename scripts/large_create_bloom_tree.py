# Usage python large_create_bloom_tree.py [Study Accession Number] [Run Accession Number] [Drive Type]

import sys
import subprocess
import pandas as pd

# Assign study and run accession numbers, kmer value and drive type to variables
study = sys.argv[1]
run = sys.argv[2]
kmer = str(20)
drive = int(sys.argv[3])
bc_size = str(2604034076)
count_size = str(518764546)
bf_size = str(878967)

# Create bloomTree directory in correct mapping directory
if drive == 0:
	newdir = '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/bloomTrees'
else:
	newdir = '/media/nickeener/External_Drive/'+study+'/bloomTrees'
subprocess.call(['mkdir', newdir])

'''# Calculate kmer frequencies of the large read file using ntcard and use to calculate appropriate bloom filter size
print('ntcard --kmer='+kmer+' --threads=8 --pref='+run+' '+run+'_*.fastq.gz')
subprocess.call(['./large_ntcard.sh', kmer, study, run])
if drive == 0:
	data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'+study+'_k'+kmer+'.hist', header=None, sep='\t', names=['1', '2'])
else:	
	data = pd.read_csv('/media/nickeener/External_Drive/'+study+'/'+run+'_k'+kmer+'.hist', header=None, sep='\t', names=['1', '2'])
bc_size = data['2'][1]
bc_size = str(int(bc_size+(bc_size*0.05)))
count_size = data['2'][1]-data['2'][2]
count_size = str(int(count_size+(count_size*0.05)))
if drive == 0:
	subprocess.call(['rm', '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'+study+'_k'+kmer+'.hist'])
else:
	subprocess.call(['rm', '/media/nickeener/External_Drive/'+study+'/'+run+'_k'+kmer+'.hist'])

# Calculate kmer frequencies of all the study's read files using ntcard and use to calculate appropriate bloom filter size
print('ntcard --kmer='+kmer+' --threads=8 --pref='+study+' *.fastq.gz')
subprocess.call(['./ntcard.sh', kmer, study])
#data = pd.read_csv('/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'+study+'_k'+kmer+'.hist', header=None, sep='\t', names=['1', '2'])
data = pd.read_csv('/media/nickeener/External_Drive/'+study+'/'+study+'_k'+kmer+'.hist', header=None, sep='\t', names=['1', '2'])
bf_size = data['2'][1]-data['2'][2]
#subprocess.call(['rm', '/home/nickeener/projects/drosophilaViruses/mapping/'+study+'/'+study+'_k'+kmer+'.hist'])
subprocess.call(['rm', '/media/nickeener/External_Drive/'+study+'/'+study+'_k'+kmer+'.hist'])
bf_size = (bf_size+int(bf_size*.05))/1000 # Add a small portion to make a slight overstimation and divide by 1000 to get it in K format
bf_size = str(bf_size)

print("bc_size = "+bc_size+"\ncount_size = "+count_size+"\nbf_size = "+bf_size)'''

# Call large_makebf.sh
subprocess.call(['./large_makebf.sh', study, run, bc_size, count_size, bf_size]) 
