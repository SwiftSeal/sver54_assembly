#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/stringtie.%j.out
#SBATCH -e logs/stringtie.%j.err

AGAT="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0"
TRANSDECODER="singularity exec -B /mnt/:/mnt/ docker://trinityrnaseq/transdecoder:latest"

mkdir -p results/stringtie

MERGED_BAM="results/star/merged.bam"

singularity exec -b /mnt/:/mnt/ docker://quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0 stringtie \
  ${MERGED_BAM} \
  -o results/stringtie/stringtie.gtf \
  -p 16

$AGAT agat_convert_sp_gxf2gxf.pl \
  -g results/stringtie/stringtie.gtf \
  -o results/stringtie/stringtie.gff

$AGAT agat_sp_extract_sequences.pl \
  -g results/stringtie/stringtie.gff \
  -f results/final_assembly/final_assembly.fa \
  -t exon \
  --merge \
  -o results/stringtie/stringtie.exon.fa

$TRANSDECODER TransDecoder.LongOrfs \
  -t results/stringtie/stringtie.exon.fa

$AGAT agat_sp_convert_sp_gxf2gff.pl \
  -g results/stringtie/stringtie.exon.fa.transdecoder_dir/longest_orfs.pep.gff3 \
  -o results/stringtie/stringtie.transdecoder.gff

$AGAT agat_sp_extract_sequences.pl \
  -g results/stringtie/stringtie.transdecoder.gff \
  -f results/final_assembly/final_assembly.fa \
  -p -t cds \
  -o results/stringtie/stringtie.transdecoder.pep.fa

$AGAT agat_convert_sp_GFF2BED.pl \
  --gff results/stringtie/stringtie.transdecoder.gff \
  --sub gene \
  -o results/stringtie/stringtie.transdecoder.bed