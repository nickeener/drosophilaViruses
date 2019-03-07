#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1
cd /media/nickeener/External_Drive/$1
echo cat
cat $2_*.fastq.gz > $2_combined.fastq.gz
echo gzip
gzip -d $2_combined.fastq.gz
echo jellyfish bc
jellyfish bc --canonical --mer-len=20 --size=$3 --threads=8 --output=$2.bc $2_combined.fastq
echo jellyfish count
jellyfish count --canonical --mer-len=20 --size=$4 --threads=8 --output=$2.jf --bc=$2.bc $2_combined.fastq
rm $2_combined.fastq
rm $2.bc
echo jellyfish dump
jellyfish dump --column --lower-count=2 $2.jf \
	| awk '{print $1}' \
	| ~/HowDeSBT/howdesbt makebf /dev/stdin --kmersin K=20 --bits=$5K --out=bloomTrees/$2.bf
rm $2.jf