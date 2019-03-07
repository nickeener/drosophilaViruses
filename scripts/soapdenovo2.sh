#!/bin/bash

SOAPdenovo-63mer all -p 8 -s /home/nickeener/projects/drosophilaViruses/scripts/soap.config -o /home/nickeener/projects/drosophilaViruses/mapping/SRP067567/assembly/K13 -K13 -N 10000;

SOAPdenovo-63mer all -p 8 -s /home/nickeener/projects/drosophilaViruses/scripts/soap.config -o /home/nickeener/projects/drosophilaViruses/mapping/SRP067567/assembly/K20 -K20 -N 10000;

SOAPdenovo-63mer all -p 8 -s /home/nickeener/projects/drosophilaViruses/scripts/soap.config -o /home/nickeener/projects/drosophilaViruses/mapping/SRP067567/assembly/K25 -K25 -N 10000;

SOAPdenovo-63mer all -p 8 -s /home/nickeener/projects/drosophilaViruses/scripts/soap.config -o /home/nickeener/projects/drosophilaViruses/mapping/SRP067567/assembly/K31 -K31 -N 10000
