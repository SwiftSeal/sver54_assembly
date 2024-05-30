#!/bin/bash

#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/stringtie.%j.out
#SBATCH -e logs/stringtie.%j.err

AGAT="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0"
TRANSDECODER="singularity exec -B /mnt/:/mnt/ docker://trinityrnaseq/transdecoder:latest"
STRINGTIE="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0 stringtie"

mkdir -p results/stringtie

MERGED_BAM="results/star/merged.bam"

SAMPLES=(
"0hr_Pinf_infection_Rep1_S11"
"0hr_Pinf_infection_Rep2_S19"
"0hr_Pinf_infection_Rep3_S27"
"24hr_Pinf_infection_Rep1_S12"
"24hr_Pinf_infection_Rep2_S20"
"24hr_Pinf_infection_Rep3_S28"
"Leaf_Rep1_S8"
"Leaf_Rep2_S16"
"Leaf_Rep3_S24"
"Root_Rep1_S10"
"Root_Rep2_S18"
"Root_Rep3_S26"
"Salt_100_mM_15_days_S34"
"Salt_100_mM_30_days_S37"
"Salt_200_mM_15_days_S35"
"Salt_control_0mM_15_days_S33"
"Salt_control_0mM_30_days_S36"
"Shoot_Rep1_S9"
"Shoot_Rep2_S17"
"Shoot_Rep3_S25"
"Temperature_stress_25C_Rep1_S14"
"Temperature_stress_25C_Rep2_S22"
"Temperature_stress_25C_Rep3_S30"
"Temperature_stress_35C_Rep1_S15"
"Temperature_stress_35C_Rep2_S23"
"Temperature_stress_35C_Rep3_S31"
"Temperature_stress_4C_Rep1_S13"
"Temperature_stress_4C_Rep2_S21"
"Temperature_stress_4C_Rep3_S29"
"Tuber_S32"
)

for SAMPLE in "${SAMPLES[@]}"; do
  $STRINGTIE \
    -p 16 \
    -o results/stringtie/${SAMPLE}.gtf\
    results/star/${SAMPLE}_Aligned.sortedByCoord.out.bam
done

$STRINGTIE \
  --merge \
  -o results/stringtie/stringtie.gtf \
  results/stringtie/*.gtf

cd results/stringtie
$TRANSDECODER $APPS/TransDecoder/util/gtf_genome_to_cdna_fasta.pl stringtie.gtf ../final_assembly/final_assembly.fa > transcripts.fa
$TRANSDECODER $APPS/TransDecoder/util/gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3
$TRANSDECODER TransDecoder.LongOrfs -t transcripts.fa
$TRANSDECODER TransDecoder.Predict -t transcripts.fa
$TRANSDECODER $APPS/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl \
  transcripts.fa.transdecoder.gff3 \
  transcripts.gff3 \
  transcripts.fa \
  > transdecoder.gff
cd ../..

$AGAT agat_convert_sp_gxf2gxf.pl \
  -g results/stringtie/transdecoder.gff \
  -o results/stringtie/final.gff

$AGAT agat_sp_extract_sequences.pl \
  -g results/stringtie/final.gff \
  -f results/final_assembly/final_assembly.fa \
  -p -t cds \
  -o results/stringtie/final.pep.fa

$AGAT agat_sp_keep_longest_isoform.pl \
    --gff results/stringtie/final.gff \
    -o results/stringtie/final.longest.gff

$AGAT agat_sp_extract_sequences.pl \
    -g results/stringtie/final.longest.gff \
    -f results/final_assembly/final_assembly.fa \
    -p -t cds \
    -o results/stringtie/final.longest.pep.fa