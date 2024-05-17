#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/braker3.%j.out
#SBATCH -e logs/braker3.%j.err

ASSEMBLY_FILE="results/final_assembly/final_assembly.fa.new.masked"
MERGED_BAM="results/star/merged.bam"

singularity exec -B /mnt/:/mnt/ $APPS/braker/braker3.sif braker.pl \
  --AUGUSTUS_CONFIG_PATH=$APPS/braker/Augustus/config/ \
  --genome=${ASSEMBLY_FILE} \
  --bam=${MERGED_BAM} \
  --prot_seq=$APPS/braker/Viridiplantae.fa \
  --softmasking \
  --threads 16 \
  --workingdir=results/braker
