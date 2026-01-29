#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/yahs.%j.out
#SBATCH -e logs/yahs.%j.err

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)
PAIRTOOLS_BAM="results/scaffolding/mapped_pairtools.bam"

$APPS/yahs/yahs \
  --no-contig-ec \
  $ASSEMBLY_PATH \
  $PAIRTOOLS_BAM \
  -o results/scaffolding/yahs.out

$APPS/yahs/juicer pre \
  -a \
  -o results/scaffolding/out_JBAT \
  results/scaffolding/yahs.out.bin \
  results/scaffolding/yahs.out_scaffolds_final.agp \
  $ASSEMBLY_PATH.fai > results/scaffolding/out_JBAT.log 2>&1

java -jar -Xmx32G $APPS/juicer_tools.1.9.9_jcuda.0.8.jar pre \
  results/scaffolding/out_JBAT.txt \
  results/scaffolding/out_JBAC.hic \
  <(cat results/scaffolding/out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')
