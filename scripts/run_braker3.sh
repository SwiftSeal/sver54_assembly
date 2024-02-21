#!/bin/bash

#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/braker3.%j.out
#SBATCH -e logs/braker3.%j.err

singularity exec -H $PWD $APPS/braker3.sif braker.pl \
        --AUGUSTUS_CONFIG_PATH=$PWD/config/config \
        --genome=results/final_assembly/solanum_verrucosum.fasta \
        --bam= \
        --prot_seq= \
        --softmasking \
        --threads {threads} \
        --workingdir=results/braker