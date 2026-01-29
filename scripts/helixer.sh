#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --gpus=1
#SBATCH --export=ALL
#SBATCH -o logs/helixer.%j.log
#SBATCH -e logs/helixer.%j.log

mkdir -p results/helixer

singularity exec --nv -B /mnt/:/mnt/ containers/helixer.sif \
  Helixer.py \
  --species solanum_verrucosum \
  --fasta-path results/assembly/cpc54.assembly.fa \
  --gff-output-path results/helixer/cpc54.helixer.gff \
  --lineage land_plant