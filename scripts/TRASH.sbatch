#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH --export=ALL
#SBATCH -o logs/trash.%j.out
#SBATCH -e logs/trash.%j.lout

mkdir -p results/TRASH

source activate TRASH
$APPS/TRASH/TRASH_run.sh results/assembly/cpc54.assembly.fa $(realpath results/TRASH) --win 10000 --m 9000 --par $SLURM_CPUS_PER_TASK

awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1, $4-1, $5}' results/TRASH/TRASH_cpc54.assembly.fa.gff > results/TRASH/TRASH_cpc54.assembly.fa.bed
