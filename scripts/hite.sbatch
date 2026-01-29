#!/bin/bash

#SBATCH -p medium
#SBATCH -c 32
#SBATCH --mem 64gb
#SBATCH --export all
#SBATCH -o logs/hite.%j.out
#SBATCH -e logs/hite.%j.out

singularity run -B /mnt/:/mnt/ docker://kanghu/hite:3.3.3 python /HiTE/main.py \
  --genome results/assembly/cpc54.assembly.fa \
  --annotate 1 \
  --thread $SLURM_CPUS_PER_TASK \
  --out_dir results/hite
