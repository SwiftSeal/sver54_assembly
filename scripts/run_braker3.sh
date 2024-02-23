#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --export=ALL
#SBATCH -o logs/braker3.%j.out
#SBATCH -e logs/braker3.%j.err

ASSEMBLY_FILE="../scratch/solanum_verrucosum/results/genome/solanum_verrucosum.fa"
MERGED_BAM="../scratch/solanum_verrucosum/results/STAR_align/Temperature_stress_4C_Rep3_Aligned.sortedByCoord.out.bam"

singularity exec -B /mnt/:/mnt/ $APPS/braker/braker3.sif braker.pl \
    --AUGUSTUS_CONFIG_PATH=$APPS/braker/Augustus/config/ \
    --genome=${ASSEMBLY_FILE} \
    --bam=${MERGED_BAM} \
    --prot_seq=$APPS/braker/Viridiplantae.fa \
    --softmasking \
    --threads 8 \
    --workingdir=results/braker
