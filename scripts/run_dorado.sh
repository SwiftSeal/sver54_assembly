#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --gpus=2
#SBATCH --export=ALL
#SBATCH -o logs/dorado.%j.out
#SBATCH -e logs/dorado.%j.err

$APPS/dorado-0.3.4-linux-x64/bin/dorado basecaller -r --emit-fastq \
$APPS/dorado-0.3.4-linux-x64/models/dna_r9.4.1_e8_sup@v3.6 \
$PROJECTS/jhi/potato/202205_Sver-ONT_Moray > assembly/dorado.fastq
