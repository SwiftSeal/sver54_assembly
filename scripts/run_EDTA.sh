#!/bin/bash

#SBATCH -p medium
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/edta.%j.out
#SBATCH -e logs/edta.%j.err

ASSEMBLY_FASTA="results/final_assembly/final_assembly.fa"

mkdir -p results/edta
cd $TMPDIR

source activate EDTA2

$APPS/EDTA/EDTA.pl --genome "$ASSEMBLY_FASTA" \
    --threads 32 \
    --anno 1

cd $SLURM_SUBMIT_DIR

cp -r $TMPDIR/* results/edta/
