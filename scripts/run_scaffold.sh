#!/bin/bash

#SBATCH -p long
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --export=ALL
#SBATCH -o logs/scaffold.%j.out
#SBATCH -e logs/scaffold.%j.err

source activate nextflow

nextflow run nextflow/scaffold.nf \
  -with-report reports/scaffold.html \
  --assembly results/quickmerge/merged_quickmerge.fasta \
  --R1_1 "/mnt/shared/projects/jhi/potato/202306_Sver-HiC/sver_54_1_S1_R1_001.fastq.gz" \
  --R1_2 "/mnt/shared/projects/jhi/potato/202306_Sver-HiC/sver_54_2_S2_R1_001.fastq.gz" \
  --R2_1 "/mnt/shared/projects/jhi/potato/202306_Sver-HiC/sver_54_1_S1_R2_001.fastq.gz" \
  --R2_2 "/mnt/shared/projects/jhi/potato/202306_Sver-HiC/sver_54_2_S2_R2_001.fastq.gz"
