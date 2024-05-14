#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/resistify_helixer.%j.out
#SBATCH -e logs/resistify_helixer.%j.err

mkdir -p results/resistify

source activate resistify

resistify --threads 8 results/helixer/helixer.pep.fa results/resistify/helixer