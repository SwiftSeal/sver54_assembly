#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/eggnog.%j.out
#SBATCH -e logs/egnnog.%j.err

source activate eggnog

emapper.py -i results/final_annotation/final_annotation.pep.fa -o results/eggnog/eggnog --cpu 16 --override