#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/windows.%j.out
#SBATCH -e logs/windows.%j.err

BASE_DIR="/mnt/shared/scratch/msmith/solanum_verrucosum"
ASSEMBLY_FASTA="$BASE_DIR/results/genome/solanum_verrucosum.fa"
GENE_BED="$BASE_DIR/results/genes/solanum_verrucosum.bed"
EDTA_GFF="$BASE_DIR/results/edta/solanum_verrucosum.fa.mod.EDTA.TEanno.gff3"

source activate windows

# get chromosome sizes file
samtools faidx
cut -f1,2 "$ASSEMBLY_FASTA.fai" > "$TMPDIR/chrom.sizes"

# make 1mb windows
bedtools makewindows \
    -g "$TMPDIR/chrom.sizes" \
    -w 1000000 \
    > "$TMPDIR/windows_1mb.bed"

# calculate methylation coverage
multiBigwigSummary BED-file \
    -b "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CG.bw" \
    --bed "$TMPDIR/windows_1mb.bed" \
    -out results/windows/methylation.npz \
    --outRaw results/windows/methylation.tab

# calculate gc content
bedtools nuc \
    -fi "$ASSEMBLY_FASTA" \
    -bed "$TMPDIR/windows_1mb.bed" \
    > results/windows/gc.bed

# calculate gene coverage
bedtools coverage \
    -a "$TMPDIR/windows_1mb.bed" \
    -b "$GENE_BED" \
    > results/windows/genes.bed

# calculate TE coverage

declare -A features=(
    ["copia"]="Copia_LTR_retrotransposon"
    ["gypsy"]="Gypsy_LTR_retrotransposon"
    ["helitron"]="helitron"
    ["tir"]="hAT_TIR_transposon Mutator_TIR_transposon PIF_Harbinger_TIR_transposon Tc1_Mariner_TIR_transposon CACTA_TIR_transposon"
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
    done < "$EDTA_GFF"
done

bedtools coverage \
    -a $TMPDIR/windows_1mb.bed \
    -b $TMPDIR/copia.bed \
    > results/windows/copia.bed

bedtools coverage \
    -a $TMPDIR/windows_1mb.bed \
    -b $TMPDIR/gypsy.bed \
    > results/windows/gypsy.bed

bedtools coverage \
    -a $TMPDIR/windows_1mb.bed \
    -b $TMPDIR/helitron.bed \
    > results/windows/helitron.bed

bedtools coverage \
    -a $TMPDIR/windows_1mb.bed \
    -b $TMPDIR/tir.bed \
    > results/windows/tir.bed
