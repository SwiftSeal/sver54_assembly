#!/bin/bash

#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/remapping.%j.out
#SBATCH -e logs/remapping.%j.err

# why write nice code when copy and paste exists

BOWTIE2="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1 bowtie2"
BOWTIE2_BUILD="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1 bowtie2-build"
SAMTOOLS="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 samtools"
VARSCAN="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/varscan:2.4.6--hdfd78af_0 varscan"

mkdir -p results/remapping

GENSEQ_PARENT_RESISTANT="-1 data/ERR2222736_1.fastq.gz -2 data/ERR2222736_2.fastq.gz"
GENSEQ_PARENT_SUSCEPTIBLE="-1 data/ERR2222737_1.fastq.gz -2 data/ERR2222737_2.fastq.gz"
GENSEQ_BULK_RESISTANT="-1 data/ERR2222738_1.fastq.gz -2 data/ERR2222738_2.fastq.gz"
GENSEQ_BULK_SUSCEPTIBLE="-1 data/ERR2222739_1.fastq.gz -2 data/ERR2222739_2.fastq.gz"

RENSEQ_PARENT_RESISTANT="-1 data/ERR2222741_1.fastq.gz -2 data/ERR2222741_2.fastq.gz"
RENSEQ_PARENT_SUSCEPTIBLE="-1 data/ERR2222742_1.fastq.gz -2 data/ERR2222742_2.fastq.gz"
RENSEQ_BULK_RESISTANT="-1 data/ERR2222744_1.fastq.gz -2 data/ERR2222744_2.fastq.gz"
RENSEQ_BULK_SUSCEPTIBLE="-1 data/ERR2222743_1.fastq.gz -2 data/ERR2222743_2.fastq.gz"

$BOWTIE2_BUILD results/final_assembly/final_assembly.fa $TMPDIR/final_assembly

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $GENSEQ_PARENT_RESISTANT \
  > results/remapping/genseq_parent_resistant.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $GENSEQ_PARENT_SUSCEPTIBLE \
  > results/remapping/genseq_parent_susceptible.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $GENSEQ_BULK_RESISTANT \
  > results/remapping/genseq_bulk_resistant.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $GENSEQ_BULK_SUSCEPTIBLE \
  > results/remapping/genseq_bulk_susceptible.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $RENSEQ_PARENT_RESISTANT \
  > results/remapping/renseq_parent_resistant.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $RENSEQ_PARENT_SUSCEPTIBLE \
  > results/remapping/renseq_parent_susceptible.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $RENSEQ_BULK_RESISTANT \
  > results/remapping/renseq_bulk_resistant.sam

bowtie2 \
  -p 4 \
  --score-min L,-0.18,-0.18 \
  --phred33 \
  --fr \
  --maxins 1000 \
  --very-sensitive \
  --no-unal \
  --no-discordant \
  -x $TMPDIR/final_assembly \
  $RENSEQ_BULK_SUSCEPTIBLE \
  > results/remapping/renseq_bulk_susceptible.sam

samtools sort -@ 4 -o results/remapping/genseq_parent_resistant.bam results/remapping/genseq_parent_resistant.sam
samtools sort -@ 4 -o results/remapping/genseq_parent_susceptible.bam results/remapping/genseq_parent_susceptible.sam
samtools sort -@ 4 -o results/remapping/genseq_bulk_resistant.bam results/remapping/genseq_bulk_resistant.sam
samtools sort -@ 4 -o results/remapping/genseq_bulk_susceptible.bam results/remapping/genseq_bulk_susceptible.sam
samtools sort -@ 4 -o results/remapping/renseq_parent_resistant.bam results/remapping/renseq_parent_resistant.sam
samtools sort -@ 4 -o results/remapping/renseq_parent_susceptible.bam results/remapping/renseq_parent_susceptible.sam
samtools sort -@ 4 -o results/remapping/renseq_bulk_resistant.bam results/remapping/renseq_bulk_resistant.sam
samtools sort -@ 4 -o results/remapping/renseq_bulk_susceptible.bam results/remapping/renseq_bulk_susceptible.sam

samtools mpileup -f results/final_assembly/final_assembly.fa results/remapping/genseq_parent_resistant.bam results/remapping/genseq_parent_susceptible.bam > results/remapping/genseq_parent.pileup
samtools mpileup -f results/final_assembly/final_assembly.fa results/remapping/genseq_bulk_resistant.bam results/remapping/genseq_bulk_susceptible.bam > results/remapping/genseq_bulk.pileup
samtools mpileup -f results/final_assembly/final_assembly.fa results/remapping/renseq_parent_resistant.bam results/remapping/renseq_parent_susceptible.bam > results/remapping/renseq_parent.pileup
samtools mpileup -f results/final_assembly/final_assembly.fa results/remapping/renseq_bulk_resistant.bam results/remapping/renseq_bulk_susceptible.bam > results/remapping/renseq_bulk.pileup

varscan mpileup2snp results/remapping/genseq_parent.pileup --output-vcf 1 --strand-filter 0 > results/remapping/genseq_parent_snp.vcf
varscan mpileup2snp results/remapping/genseq_bulk.pileup --output-vcf 1 --strand-filter 0 > results/remapping/genseq_bulk_snp.vcf
varscan mpileup2snp results/remapping/renseq_parent.pileup --output-vcf 1 --strand-filter 0 > results/remapping/renseq_parent_snp.vcf
varscan mpileup2snp results/remapping/renseq_bulk.pileup --output-vcf 1 --strand-filter 0 > results/remapping/renseq_bulk_snp.vcf