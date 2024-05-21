#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/merge.%j.out
#SBATCH -e logs/merge.%j.err

source activate sequtils

mkdir -p results/final_annotation

# Complement BRAKER3 annotation with Helixer NLRs

awk '{print $1}' results/resistify/helixer/results.tsv \
   | tail -n +1 \
   | sort \
   | uniq \
   > $TMPDIR/helixer_resistify.genes

agat_sp_filter_feature_from_keep_list.pl \
  --gff results/helixer/helixer.gff \
  --keep_list $TMPDIR/helixer_resistify.genes \
  -o $TMPDIR/helixer_nlrs.gff

agat_sp_complement_annotations.pl \
  --ref results/braker/braker.gff \
  --add $TMPDIR/helixer_nlrs.gff \
  -o results/final_annotation/final_annotation.gff

# Extract sequences

agat_sp_extract_sequences.pl \
  -g results/final_annotation/final_annotation.gff \
  -f results/final_assembly/final_assembly.fa \
  -p -t cds \
  -o results/final_annotation/final_annotation.pep.fa

agat_sp_extract_sequences.pl \
  -g results/final_annotation/final_annotation.gff \
  -f results/final_assembly/final_assembly.fa \
  -t cds \
  -o results/final_annotation/final_annotation.cds.fa

agat_convert_sp_gff2bed \
  --gff results/final_annotation/final_annotation.gff \
  --sub gene \
  -o results/final_annotation/final_annotation.bed

agat_sp_keep_longest_isoform.pl \
  --gff results/final_annotation/final_annotation.gff \
  -o results/final_annotation/final_annotation.longest.gff

agat_sp_extract_sequences.pl \
  -g results/final_annotation/final_annotation.longest.gff \
  -f results/final_assembly/final_assembly.fa \
  -p -t cds \
  -o results/final_annotation/final_annotation.longest.pep.fa

agat_sp_extract_sequences.pl \
  -g results/final_annotation/final_annotation.longest.gff \
  -f results/final_assembly/final_assembly.fa \
  -t cds \
  -o results/final_annotation/final_annotation.longest.cds.fa

agat_convert_sp_gff2bed.pl \
  --gff results/final_annotation/final_annotation.longest.gff \
  --sub gene \
  -o results/final_annotation/final_annotation.longest.bed