#!/bin/bash

# Convert fasta files to Bloom filter bit vectors
cd ~/projects/drosophilaViruses/mapping/$1
cd /media/nickeener/External_Drive/$1
cat $2*.fastq.gz > $2_combined.fastq.gz
gzip -dc $2_combined.fastq.gz > bloomTrees/$2_combined.fastq
rm $2_combined.fastq.gz
cd bloomTrees
~/tools/HowDeSBT/howdesbt makebf --k=$3 --min=2 --bits=$4K --threads=8 $2_combined.fastq --out=$2.bf
rm $2_combined.fastq