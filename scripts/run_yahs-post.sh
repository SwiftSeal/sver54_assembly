#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/yahs-post.%j.out
#SBATCH -e logs/yahs-post.%j.err

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)

$APPS/yahs/juicer post \
  -o results/scaffolding/out_JBAT \
  results/scaffolding/out_JBAT.review.assembly \
  results/scaffolding/out_JBAT.liftover.agp \
  $ASSEMBLY_PATH
