#!/bin/bash

#SBATCH -p himem
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/trash.%j.out
#SBATCH -e logs/trash.%j.err

mkdir -p results/TRASH

source activate TRASH
$APPS/TRASH/TRASH_run.sh results/final_assembly/final_assembly.fa $(realpath results/TRASH) --win 10000 --m 9000 --par 8