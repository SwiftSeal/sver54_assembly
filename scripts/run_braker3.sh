#!/bin/bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/braker3.%j.out
#SBATCH -e logs/braker3.%j.err

ASSEMBLY_FILE="results/final_assembly/final_assembly.fa.new.masked"
MERGED_BAM="results/star/merged.bam"
AGAT_SP_EXTRACT_SEQUENCES="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0 agat_sp_extract_sequences.pl"
AGAT_CONVERT_SP_GXF2GXF="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0 agat_convert_sp_gxf2gxf.pl"
AGAT_CONVERT_SP_GFF2BED="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0 agat_convert_sp_gff2bed.pl"

singularity exec -B /mnt/:/mnt/ $APPS/braker/braker3.sif braker.pl \
  --AUGUSTUS_CONFIG_PATH=$APPS/braker/Augustus/config/ \
  --genome=${ASSEMBLY_FILE} \
  --bam=${MERGED_BAM} \
  --prot_seq=$APPS/braker/Viridiplantae.fa \
  --softmasking \
  --threads 16 \
  --workingdir=results/braker

$AGAT_CONVERT_SP_GXF2GXF -g results/braker/braker.gtf -o results/braker/braker.gff
$AGAT_SP_EXTRACT_SEQUENCES -g results/braker/braker.gff -f results/final_assembly/final_assembly.fa -p -t cds -o results/braker/braker.pep.fa
$AGAT_SP_EXTRACT_SEQUENCES -g results/braker/braker.gff -f results/final_assembly/final_assembly.fa -t cds -o results/braker/braker.cds.fa
$AGAT_CONVERT_SP_GFF2BED --gff results/braker/braker.gff --sub gene -o results/braker/braker.bed
