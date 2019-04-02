#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1
cd /media/nickeener/External_Drive/$1

zcat $2_*.fastq.gz | jellyfish count /dev/fd/0 --mer-len=20 --threads=8 --size=834090000 -C --output=$2.jf
jellyfish dump $2.jf > dump.fa
rm $2.jf
