#!/bin/bash

# Only mapped reads (SAM output)
bowtie2 -t --very-sensitive-local --no-unal -p 8 -x $1 -1 $2 -2 $3 -S $4/$5.sam 2> $4/$5.log
# Only unmapped reads (fastq.gz output)
#bowtie2 -t --sensitive-local --un-conc-gz $4/$5.fastq.gz -p 8 -x $1 -1 $2 -2 $3 > /dev/null