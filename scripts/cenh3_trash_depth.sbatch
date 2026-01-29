#!/bin/env -S pixi exec -s mosdepth==0.3.12 -- bash

#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH --export=ALL
#SBATCH -o logs/cenh3_trash_depth.%j.log
#SBATCH -e logs/cenh3_trash_depth.%j.log

mosdepth -b results/TRASH/TRASH_cpc54.assembly.fa.bed results/cenh3/trash results/cenh3/SRR18548893.bam
