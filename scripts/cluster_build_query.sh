#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1/bloomTrees
cd /media/nickeener/External_Drive/$1/bloomTrees

# Create a tree topology
ls *.bf > leafnames
~/tools/HowDeSBT/howdesbt cluster --list=leafnames --bits=$2K --tree=union.sbt --nodename=node{number} --keepallnodes
rm leafnames

# Build the HowDeSBT nodes
~/tools/HowDeSBT/howdesbt build --HowDe --tree=union.sbt --outtree=howde.sbt

# Run queries
~/tools/HowDeSBT/howdesbt query --sort --time --threshold=0.01 --tree=howde.sbt  /home/nickeener/projects/drosophilaViruses/viralSeqs/twyford_full.fasta > queries_twyford_full.dat
~/tools/HowDeSBT/howdesbt query --sort --time --threshold=0.01 --tree=howde.sbt  /home/nickeener/projects/drosophilaViruses/viralSeqs/twyford.fasta > queries_twyford.dat
#~/HowDeSBT/howdesbt query --sort --time --threshold=0.05 --tree=howde.sbt  /home/nickeener/projects/drosophilaViruses/viralSeqs/Emuscae_TSA.fasta > queries_Emuscae.dat