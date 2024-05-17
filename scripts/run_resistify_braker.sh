#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/resistify_braker.%j.out
#SBATCH -e logs/resistify_braker.%j.err

mkdir -p results/resistify

source activate resistify

resistify --threads 8 results/braker/braker.pep.fa results/resistify/braker