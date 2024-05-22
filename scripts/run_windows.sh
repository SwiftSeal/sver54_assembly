#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/windows.%j.out
#SBATCH -e logs/windows.%j.err

MULTIBIGWIGSUMMARY="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 multiBigwigSummary"
SAMTOOLS="singularity exec -B /mnt/:/mnt/ https://depot.galaxyproject.org/singularity/samtools:1.19.1--h50ea8bc_0 samtools"
BEDTOOLS="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1 bedtools"

ASSEMBLY_FASTA="results/final_assembly/final_assembly.fa"
GENE_BED="results/final_annotation/final_annotation.longest.bed"
EARLGREY_GFF="results/earlgrey/solanum_verrucosum_EarlGrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.filteredRepeats.gff"

# get chromosome sizes file
$SAMTOOLS faidx "$ASSEMBLY_FASTA"
cut -f1,2 "$ASSEMBLY_FASTA.fai" > "$TMPDIR/chrom.sizes"

# make 1mb windows
$BEDTOOLS makewindows \
  -g "$TMPDIR/chrom.sizes" \
  -w 1000000 \
  > "$TMPDIR/windows_1mb.bed"

# calculate methylation coverage
#$MULTIBIGWIGSUMMARY BED-file \
#  -b \
#    "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CG.bw" \
#    "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CHG.bw" \
#    "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CHH.bw" \
#  --BED "$TMPDIR/windows_1mb.bed" \
#  -out results/windows/methylation.npz \
#  --outRaw results/windows/methylation.tab

# calculate gc content
$BEDTOOLS nuc \
  -fi "$ASSEMBLY_FASTA" \
  -bed "$TMPDIR/windows_1mb.bed" \
  > results/windows/gc.bed

# calculate gene coverage
$BEDTOOLS coverage \
  -a "$TMPDIR/windows_1mb.bed" \
  -b "$GENE_BED" \
  > results/windows/genes.bed

# calculate TE coverage

declare -A features=(
  ["Ty1"]="LTR/Copia"
  ["Ty3"]="LTR/Gypse"
  ["Helitron"]="RC/Helitron"
)

for feature in "${!features[@]}"; do
  while IFS=$'\t' read -r col1 col2 col3 col4 col5; do
    if [[ $col1 == "#" ]]; then
      continue
    fi
    for f in ${features[$feature]}; do
      if [[ $col3 == $f ]]; then
        echo -e "$col1\t$col4\t$col5" >> "$TMPDIR/${feature}.bed"
      fi
    done
  done < "$EARLGREY_GFF"
done

$BEDTOOLS coverage \
  -a $TMPDIR/windows_1mb.bed \
  -b $TMPDIR/Ty1.bed \
  > results/windows/Ty1.bed

$BEDTOOLS coverage \
  -a $TMPDIR/windows_1mb.bed \
  -b $TMPDIR/Ty3.bed \
  > results/windows/Ty3.bed

$BEDTOOLS coverage \
  -a $TMPDIR/windows_1mb.bed \
  -b $TMPDIR/Helitron.bed \
  > results/windows/Helitron.bed
