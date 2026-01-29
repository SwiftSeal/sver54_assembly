#!/bin/bash

#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/edta.%j.log
#SBATCH -e logs/edta.%j.log

ASSEMBLY_FASTA="$(realpath results/assembly/cpc54.assembly.fa)"

# EDTA runs only from local environment
mkdir -p results/edta
cd results/edta

singularity exec -B /mnt/:/mnt/ ../../containers/edta.sif \
    EDTA.pl \
    --genome $ASSEMBLY_FASTA \
    --sensitive 1 \
    --threads 32 \
    --anno 1
