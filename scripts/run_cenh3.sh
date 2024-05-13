#!/bin/bash

#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/cenh3.%j.out
#SBATCH -e logs/cenh3.%j.err

MINIMAP2="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/minimap2:2.28--he4a0461_1 minimap2"
SAMTOOLS="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 samtools"
BAMCOVERAGE="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bamCoverage"

mkdir -p results/cenh3

$MINIMAP2 -t 4 -ax sr results/final_assembly/final_assembly.fa data/SRR18548893.fastq.gz > results/cenh3/aln.sam

$SAMTOOLS view -b -o results/cenh3/aln.bam results/cenh3/aln.sam
$BAMCOVERAGE -b results/cenh3/aln.bam -o results/cenh3/aln.bw
