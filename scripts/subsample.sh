#!/bin/bash

cd ~/projects/drosophilaViruses/mapping/$1
#cd /media/nickeener/External_Drive/$1

/home/nickeener/bbmap/reformat.sh in=$2 in2=$3 out=subsample$4/$2 out2=subsample$4/$3 samplerate=$4