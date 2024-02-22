#!/bin/bash

#SBATCH -p short
#SBATCH -c 64
#SBATCH --mem=128G
#SBATCH --export=ALL
#SBATCH -o logs/edta.%j.out
#SBATCH -e logs/edta.%j.err

ASSEMBLY_FASTA="/mnt/shared/scratch/msmith/solanum_verrucosum/results/genome/solanum_verrucosum.fa"

mkdir -p results/edta
cd $TMPDIR

source activate EDTA2

$APPS/EDTA/EDTA.pl --genome "$ASSEMBLY_FASTA" \
    --threads 64 \
    --anno 1

cd -

cp -r $TMPDIR/* results/edta/
