#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/tesorter.%j.out
#SBATCH -e logs/tesorter.%j.err

source activate tesorter

TEsorter -p 8 ~/scratch/solanum_verrucosum/results/edta/solanum_verrucosum.fa.mod.EDTA.final/solanum_verrucosum.fa.mod.EDTA.intact.fa
