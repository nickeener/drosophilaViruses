#!/bin/bash

cd /home/nickeener/projects/drosophilaViruses

mkdir mapping/SRP133370/assembly
mkdir mapping/SRP133370/assembly/K$1

SOAPdenovo-63mer all -p 1 -s scripts/soap.config -o mapping/SRP133370/assembly/K$1 -K $1 -N 10000;
