#!/bin/bash

zgrep -A -B 1 '<sequence>' <fastq file> | sed '/^--$/d' > out.fq