#!/bin/bash

cd /home/nickeener/projects/drosophilaViruses/mapping/$2
cd /media/nickeener/External_Drive/$2
ntcard --kmer=$1 --threads=8 --pref=$2 *.fastq.gz
