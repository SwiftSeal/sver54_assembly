#!/bin/bash

#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/compleasm.%j.out
#SBATCH -e logs/compleasm.%j.err

source activate compleasm

compleasm run \
  -a results/final_assembly/final_assembly.fa \
  -o results/compleasm \
  -t 4 \
  -l viridiplantae_odb10