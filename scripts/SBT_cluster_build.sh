#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1/bloomTrees
cd /media/nickeener/External_Drive/$1/bloomTrees

# Create a tree topology
ls *.bf > leafnames
~/tools/HowDeSBT/howdesbt cluster --list=leafnames --bits=$2K --tree=union.sbt --nodename=node{number} --keepallnodes
rm leafnames

# Build the HowDeSBT nodes
~/tools/HowDeSBT/howdesbt build --HowDe --tree=union.sbt --outtree=howde.sbt
