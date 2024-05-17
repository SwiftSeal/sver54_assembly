#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/tesorter.%j.out
#SBATCH -e logs/tesorter.%j.err

mkdir -p results/tesorter

source activate tesorter

TEsorter -p 8 results/edta/final_assembly.fa.mod.EDTA.TElib.fa -db rexdb-plant -pre results/tesorter/tesorter
