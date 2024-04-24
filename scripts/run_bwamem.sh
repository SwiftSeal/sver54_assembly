#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/pairtools.%j.out
#SBATCH -e logs/pairtools.%j.err

BWA="singularity exec -B /mnt/:/mnt/ https://depot.galaxyproject.org/singularity/bwa:0.7.16--pl5.22.0_0 bwa"

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)
READ_PATH="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/hic"
R1_1="$READ_PATH/sver_54_1_S1_R1_001.fastq.gz"
R1_2="$READ_PATH/sver_54_2_S2_R1_001.fastq.gz"
R2_1="$READ_PATH/sver_54_1_S1_R2_001.fastq.gz"
R2_2="$READ_PATH/sver_54_2_S2_R2_001.fastq.gz"

mkdir -p results/scaffolding
cd $TMPDIR

$BWA index $ASSEMBLY_PATH

$BWA mem \
  -5SP \
  -T0 \
  -t 16 \
  $ASSEMBLY_PATH \
  <(zcat $R1_1 $R1_2) <(zcat $R2_1 $R2_2) \
  -o aligned.sam