#!/bin/bash

#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/resistify_final.%j.out
#SBATCH -e logs/resistify_final.%j.err

mkdir -p results/resistify

source activate resistify

resistify --threads 16 results/final_annotation/final_annotation.longest.pep.fa results/resistify/final