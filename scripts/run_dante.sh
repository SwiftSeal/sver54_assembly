#!/bin/bash

#SBATCH -p short
#SBATCH -c 32
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/dante.%j.out
#SBATCH -e logs/dante.%j.err

mkdir -p results/dante

source activate dante_ltr

dante \
  -q results/final_assembly/final_assembly.fa \
  -o results/dante/dante_output.gff3 \
  -c 32

dante_ltr \
  -g results/dante/dante_output.gff3 \
  -s results/final_assembly/final_assembly.fa \
  -o results/dante/dante_ltr_annotation \
  -c 32 \
  -M 1

dante_ltr_to_library \
  -g results/dante/dante_ltr_annotation.gff3 \
  -s results/final_assembly/final_assembly.fa \
  -o results/dante/dante_ltr_library.fa \
  -c 32