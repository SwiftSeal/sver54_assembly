#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/STAR.%j.out
#SBATCH -e logs/STAR.%j.err

ASSEMBLY_FASTA="results/genome/solanum_verrucosum.fa"
RNA_SEQ_DIR="/mnt/shared/projects/jhi/potato/202212_Sver-RNAseq/"
R1_SUFFIX="_R1_001.fastq.gz"
R2_SUFFIX="_R2_001.fastq.gz"
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

source activate star

STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $TMPDIR/index \
    --genomeFastaFiles $ASSEMBLY_FASTA

for SAMPLE in "${SAMPLES[@]}"; do
    R1="$RNA_SEQ_DIR/${SAMPLE}${R1_SUFFIX}"
    R2="$RNA_SEQ_DIR/${SAMPLE}${R2_SUFFIX}"
    STAR \
    --runThreadN 16 \
    --genomeDir $TMPDIR/index \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --readFilesIn $R1 $R2 \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$TMPDIR/${SAMPLE}_" \
    --alignIntronMin 60 \
    --alignIntronMax 15000 \
    --alignMatesGapMax 2000 \
    --alignEndsType Local \
    --alignSoftClipAtReferenceEnds No \
    --outSAMprimaryFlag AllBestScore \
    --outFilterMismatchNoverLmax 0.02 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 1 \
    --outFilterMatchNmin 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMultimapNmax 15 \
    --alignTranscriptsPerReadNmax 30000 \
    --alignSJoverhangMin 7 \
    --alignSJDBoverhangMin 7 \
    --alignSJstitchMismatchNmax 0 1 0 0
done

mkdir -p results/star

samtools merge \
    -@ 16 \
    -o results/star/merged.bam \
    $TMPDIR/*_Aligned.sortedByCoord.out.bam

mv $TMPDIR/*_Aligned.sortedByCoord.out.bam results/star/
