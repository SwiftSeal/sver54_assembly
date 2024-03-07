#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --gpus=2
#SBATCH --export=ALL
#SBATCH -o logs/guppy.%j.out
#SBATCH -e logs/guppy.%j.err

$APPS/ont-guppy/bin/guppy_basecaller \
  -i $PROJECTS/jhi/potato/202205_Sver-ONT_Moray \
  -s results/guppy/ \
  -c dna_r9.4.1_450bps_sup_plant.cfg \
  --device "auto" \
  --recursive

cat results/guppy/pass/*.fastq > results/guppy/reads.fastq
gzip results/guppy/reads.fastq
