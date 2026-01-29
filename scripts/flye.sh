#!/bin/bash

#SBATCH -p himem
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --export=ALL
#SBATCH -o logs/flye.%j.log
#SBATCH -e logs/flye.%j.log

cat results/basecalled/reads_*.fq.gz > $TMPDIR/reads.fq.gz

singularity exec -B /mnt/:/mnt containers/flye.sif \
    flye --nano-hq $TMPDIR/reads.fq.gz --out-dir results/flye --threads $SLURM_CPUS_PER_TASK
