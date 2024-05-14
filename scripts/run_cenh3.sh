#!/bin/bash

#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/cenh3.%j.out
#SBATCH -e logs/cenh3.%j.err

BOWTIE2_BUILD="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1 bowtie2-build"
BOWTIE2="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1 bowtie2"
SAMTOOLS="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 samtools"
BAMCOVERAGE="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bamCoverage"

mkdir -p results/cenh3

$BOWTIE2_BUILD results/final_assembly/final_assembly.fa $TMPDIR/final_assembly
$BOWTIE2 -p 4 -x $TMPDIR/final_assembly -U data/SRR18548893.fastq.gz -S results/cenh3/aln.sam
$SAMTOOLS sort -o results/cenh3/aln.bam results/cenh3/aln.sam
$SAMTOOLS index results/cenh3/aln.bam
$BAMCOVERAGE -b results/cenh3/aln.bam -o results/cenh3/aln.bw
