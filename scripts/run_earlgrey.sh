#!/bin/bash

#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/earlgrey.%j.out
#SBATCH -e logs/earlgrey.%j.err

source activate earlgrey

earlGrey -g results/final_assembly/final_assembly.fa -s solanum_verrucosum -o results/earlgrey -t 32