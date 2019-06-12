#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1/bloomTrees
cd /media/nickeener/External_Drive/$1/bloomTrees

~/tools/HowDeSBT/howdesbt query --sort --time --threshold=0.01 --tree=howde.sbt  /home/nickeener/projects/drosophilaViruses/viralSeqs/$2.fasta > queries_$2.dat
