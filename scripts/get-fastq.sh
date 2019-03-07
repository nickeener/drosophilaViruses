#!/bin/bash
# Change wget command when switching between SRR and ERR runs

wget -P $3 ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$2/$1/$1.sra &&
fastq-dump -O $3 --split-files --gzip /media/nickeener/External_Drive/$4/$1.sra &&
rm /media/nickeener/External_Drive/$4/$1.sra