#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1/bloomTrees
cd /media/nickeener/External_Drive/$1/bloomTrees

/home/nickeener/Downloads/HowDeSBT/howdesbt querybf --filter=$2 $3 > $4.kmer