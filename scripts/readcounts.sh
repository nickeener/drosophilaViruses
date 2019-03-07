#!/bin/bash

cd /home/nickeener/projects/drosophilaViruses/mapping/$1
cd /media/nickeener/External_Drive/$1

zcat $2 | wc -l