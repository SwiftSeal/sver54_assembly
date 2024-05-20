#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/salmon.%j.out
#SBATCH -e logs/salmon.%j.err

ASSEMBLY_FASTA="results/final_assembly/final_assembly.fa"
CDS="results/braker/braker.cds.fa"
RNA_SEQ_DIR="/mnt/shared/projects/jhi/potato/202212_Sver-RNAseq/"
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
)
SALMON="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/salmon:1.10.3--hecfa306_0 salmon"

mkdir -p results/salmon

grep "^>" $ASSEMBLY_FASTA | cut -d " " -f 1 > $TMPDIR/decoys.txt
sed -i 's/>//g' $TMPDIR/decoys.txt
cat $CDS $ASSEMBLY_FASTA > $TMPDIR/transcripts.fa


$SALMON index -t $TMPDIR/transcripts.fa -d $TMPDIR/decoys.txt -i $TMPDIR/salmon_index -p 8

for SAMPLE in "${SAMPLES[@]}"; do
    R1="$RNA_SEQ_DIR/${SAMPLE}_R1_001.fastq.gz"
    R2="$RNA_SEQ_DIR/${SAMPLE}_R2_001.fastq.gz"
    $SALMON quant -i $TMPDIR/salmon_index -l A -1 $R1 -2 $R2 -p 8 -o results/salmon/$SAMPLE
done