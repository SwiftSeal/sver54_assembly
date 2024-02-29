#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/braker3.%j.out
#SBATCH -e logs/braker3.%j.err

ASSEMBLY_FILE="../solanum_verrucosum/results/genome/solanum_verrucosum.fa"
MERGED_BAM="../solanum_verrucosum/results/STAR_align/Temperature_stress_4C_Rep3_Aligned.sortedByCoord.out.bam"

singularity exec -B /mnt/:/mnt/ $APPS/braker/braker3.sif braker.pl \
  --AUGUSTUS_CONFIG_PATH=$APPS/braker/Augustus/config/ \
  --genome=${ASSEMBLY_FILE} \
  --bam=${MERGED_BAM} \
  --prot_seq=$APPS/braker/Viridiplantae.fa \
  --softmasking \
  --threads 16 \
  --workingdir=results/braker
