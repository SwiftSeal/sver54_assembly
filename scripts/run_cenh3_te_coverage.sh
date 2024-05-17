#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/cenh3.%j.out
#SBATCH -e logs/cenh3.%j.err

$BEDTOOLS="singularity exec -B /mnt/:/mnt/ quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1 bedtools"

$BEDTOOLS coverage \
  -a results/earlgrey/solanum_verrucosum_EarlGrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.filteredRepeats.gff \
  -b results/cenh3/aln.bam \
  -d \
  > results/cenh3/cenh3_te_coverage.bed

conda activate R

Rscript - << EOF
library(data.table)
dt <- fread("results/cenh3/cenh3_te_coverage.bed", header=FALSE)
dt[, element := sub(".*ID=([^;]+);.*", "\\1", V9)]
dt <- dt[, .(mean_coverage = mean(V11)), by = element]
fwrite(dt, "results/cenh3/cenh3_te_coverage_mean.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
EOF