#!/bin/env -S pixi exec -c conda-forge -c bioconda --spec bedtools==2.31.1 --spec samtools==1.20 --spec deeptools==3.5.5 -- bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH --export ALL
#SBATCH -o logs/windows.%j.out
#SBATCH -e logs/windows.%j.out

ASSEMBLY_FASTA="results/assembly/cpc54.assembly.fa"
GENE_BED="results/annotation/cpc54.gene.bed"
TE_GFF="results/annotation/cpc54.te.all.gff"

mkdir -p results/windows

# make 1mb windows
bedtools makewindows \
  -g $ASSEMBLY_FASTA.sizes \
  -w 1000000 \
  > results/windows/windows_1mb.bed

# calculate methylation coverage
multiBigwigSummary BED-file \
  -b \
    results/methylation/cpc54.cg.bw \
    results/methylation/cpc54.chg.bw \
    results/methylation/cpc54.chh.bw \
  --BED results/windows/windows_1mb.bed \
  -out results/windows/methylation.npz \
  --outRaw results/windows/methylation.tab

# calculate gc content
bedtools nuc \
  -fi $ASSEMBLY_FASTA \
  -bed results/windows/windows_1mb.bed \
  > results/windows/gc.bed

# calculate gene coverage
bedtools coverage \
  -a results/windows/windows_1mb.bed \
  -b $GENE_BED \
  > results/windows/genes.bed

# calculate TE coverage
declare -A features=(
  ["Ty1"]="LTR/Copia"
  ["Ty3"]="LTR/Gypsy"
  ["Helitron"]="RC/Helitron"
)

for feature in "${!features[@]}"; do
  awk -v target="${features[$feature]}" '
    BEGIN { OFS="\t" }
    /^#/ { next }                   # Skip GFF header lines
    $3 == target {                  # Exact match for Column 3
        print $1, $4-1, $5 
    }' "$TE_GFF" > "results/windows/${feature}.extracted.bed"
done

bedtools coverage \
  -a results/windows/windows_1mb.bed \
  -b results/windows/Ty1.extracted.bed \
  > results/windows/Ty1.bed

bedtools coverage \
  -a results/windows/windows_1mb.bed \
  -b results/windows/Ty3.extracted.bed \
  > results/windows/Ty3.bed

