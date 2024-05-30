#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=4G
#SBATCH --export=ALL
#SBATCH -o logs/remapping.%j.out
#SBATCH -e logs/remapping.%j.err

# why write nice code when copy and paste exists

BOWTIE2="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1 bowtie2"
BOWTIE2_BUILD="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1 bowtie2-build"
SAMTOOLS="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 samtools"
VARSCAN="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/varscan:2.4.6--hdfd78af_0 varscan"
BCFTOOLS="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0 bcftools"

mkdir -p results/remapping

GENSEQ_PARENT_RESISTANT="ERR2222736"
GENSEQ_PARENT_SUSCEPTIBLE="ERR2222737"
GENSEQ_BULK_RESISTANT="ERR2222738"
GENSEQ_BULK_SUSCEPTIBLE="ERR2222739"
RENSEQ_PARENT_RESISTANT="ERR2222741"
RENSEQ_PARENT_SUSCEPTIBLE="ERR222274"
RENSEQ_BULK_RESISTANT="ERR2222744"
RENSEQ_BULK_SUSCEPTIBLE="ERR2222743"

SAMPLES=($GENSEQ_PARENT_RESISTANT $GENSEQ_PARENT_SUSCEPTIBLE $GENSEQ_BULK_RESISTANT $GENSEQ_BULK_SUSCEPTIBLE $RENSEQ_PARENT_RESISTANT $RENSEQ_PARENT_SUSCEPTIBLE $RENSEQ_BULK_RESISTANT $RENSEQ_BULK_SUSCEPTIBLE)

$BOWTIE2_BUILD results/final_assembly/final_assembly.fa $TMPDIR/final_assembly

for SAMPLE in ${SAMPLES[@]}; do
  $BOWTIE2 -p 8 -x $TMPDIR/final_assembly -1 data/${SAMPLE}_1.fastq -2 data/${SAMPLE}_2.fastq -S results/remapping/${SAMPLE}.sam
  $SAMTOOLS sort -@ 8 -o results/remapping/${SAMPLE}.bam results/remapping/${SAMPLE}.sam
  $SAMTOOLS mpileup -f results/final_assembly/final_assembly.fa results/remapping/${SAMPLE}.bam > results/remapping/${SAMPLE}.pileup
  $VARSCAN mpileup2snp results/remapping/${SAMPLE}.pileup --output-vcf 1 --strand-filter 0 > results/remapping/${SAMPLE}.vcf
  # filter based on expected ratio for each sample
  if [[ $SAMPLE == *PARENT_RESISTANT* ]]; then
    $BCFTOOLS filter -i 'AF < 0.1' results/remapping/${SAMPLE}.vcf > results/remapping/${SAMPLE}.filtered.vcf
  elif [[ $SAMPLE == *PARENT_SUSCEPTIBLE* ]]; then
    $BCFTOOLS filter -i 'AF > 0.9' results/remapping/${SAMPLE}.vcf > results/remapping/${SAMPLE}.filtered.vcf
  elif [[ $SAMPLE == *BULK_RESISTANT* ]]; then
    $BCFTOOLS filter -i 'AF > 0.4 & AF < 0.6' results/remapping/${SAMPLE}.vcf > results/remapping/${SAMPLE}.filtered.vcf
  elif [[ $SAMPLE == *BULK_SUSCEPTIBLE* ]]; then
    $BCFTOOLS filter -i 'AF > 0.9' results/remapping/${SAMPLE}.vcf > results/remapping/${SAMPLE}.filtered.vcf
  fi
done

$BCFTOOLS filter -i 'AF > 0.9' results/remapping/${GENSEQ_PARENT_SUSCEPTIBLE}.vcf > results/remapping/${GENSEQ_PARENT_SUSCEPTIBLE}.filtered.vcf
$BCFTOOLS filter -i 'AF < 0.1' results/remapping/${GENSEQ_PARENT_RESISTANT}.vcf > results/remapping/${GENSEQ_PARENT_RESISTANT}.filtered.vcf
$BCFTOOLS filter -i 'AF > 0.9' results/remapping/${RENSEQ_PARENT_SUSCEPTIBLE}.vcf > results/remapping/${RENSEQ_PARENT_SUSCEPTIBLE}.filtered.vcf
$BCFTOOLS filter -i 'AF < 0.1' results/remapping/${RENSEQ_PARENT_RESISTANT}.vcf > results/remapping/${RENSEQ_PARENT_RESISTANT}.filtered.vcf
$BCFTOOLS filter -i 'AF > 0.4 & AF < 0.6' results/remapping/${GENSEQ_BULK_RESISTANT}.vcf > results/remapping/${GENSEQ_BULK_RESISTANT}.filtered.vcf
$BCFTOOLS filter -i 'AF > 0.9' results/remapping/${GENSEQ_BULK_SUSCEPTIBLE}.vcf > results/remapping/${GENSEQ_BULK_SUSCEPTIBLE}.filtered.vcf
$BCFTOOLS filter -i 'AF > 0.4 & AF < 0.6' results/remapping/${RENSEQ_BULK_RESISTANT}.vcf > results/remapping/${RENSEQ_BULK_RESISTANT}.filtered.vcf
$BCFTOOLS filter -i 'AF > 0.9' results/remapping/${RENSEQ_BULK_SUSCEPTIBLE}.vcf > results/remapping/${RENSEQ_BULK_SUSCEPTIBLE}.filtered.vcf

bedtools intersect -A results/remapping/${RENSEQ_PARENT_SUSCEPTIBLE}.vcf