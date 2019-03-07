#!/bin/bash

bowtie2 -t --very-sensitive --no-unal -p 8 -x $1 -1 $2 -2 $3 -S $4/$5.sam