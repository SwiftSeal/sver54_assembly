#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/tesorter.%j.out
#SBATCH -e logs/tesorter.%j.err

mkdir -p results/tesorter

source activate tesorter

TEsorter \
  -p 8 \
  results/earlgrey/solanum_verrucosum_EarlGrey/solanum_verrucosum_summaryFiles/solanum_verrucosum-families.fa.strained \
  -db rexdb-plant \
  -pre results/tesorter/earlgrey
