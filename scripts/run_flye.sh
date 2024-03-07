#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/flye.%j.out
#SBATCH -e logs/flye.%j.err

source activate flye

flye --nano-hq results/guppy/reads.fastq --out-dir results/flye --threads 16
