#!/bin/bash

# This script is for testing out cd-hit on the collection of fasta files located in /projects/microprot/data/

cd-hit -i test5.fasta -o test100 -c 1.00 -n 5
cd-hit -i test5.fasta -o test90 -c 0.9 -n 5
