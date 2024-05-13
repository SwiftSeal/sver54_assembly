#!/bin/bash

#SBATCH -p himem
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/trash.%j.out
#SBATCH -e logs/trash.%j.err

mkdir -p results/TRASH

source activate TRASH
$APPS/TRASH/TRASH_run.sh results/final_assembly/final_assembly.fa $(realpath results/TRASH) --win 5000 -m 4900 --par 8