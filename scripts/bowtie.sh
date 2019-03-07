#!/bin/bash

bowtie -t --no-unal -y -p 8 $1 -1 $2 -2 $3 -S $4/$5.sam