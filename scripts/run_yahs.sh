#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/yahs.%j.out
#SBATCH -e logs/yahs.%j.err

mkdir -p results/yahs

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)
PAIRTOOLS_BAM="results/scaffolding/mapped_pairtools.bam"

$APPS/yahs/yahs \
  --no-contig-ec \
  $ASSEMBLY_PATH \
  $PAIRTOOLS_BAM \
  -o results/yahs/yahs.out

$APPS/yahs/juicer pre \
  -a \
  -o results/yahs/out_JBAT \
  results/yahs/yahs.out.bin \
  results/yahs/yahs.out_scaffolds_final.agp \
  $ASSEMBLY_PATH.fai > results/yahs/out_JBAT.log 2>&1

java -jar -Xmx32G $APPS/juicer_tools.1.9.9_jcuda.0.8.jar pre \
  results/yahs/out_JBAT.txt \
  results/yahs/out_JBAC.hic \
  <(cat results/yahs/out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')
