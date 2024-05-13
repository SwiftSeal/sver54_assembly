#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --gpus=1
#SBATCH --export=ALL
#SBATCH -o logs/helixer.%j.out
#SBATCH -e logs/helixer.%j.err

AGAT_SP_EXTRACT_SEQUENCES="singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0 agat_sp_extract_sequences.pl"

singularity exec --nv -B /mnt/:/mnt/ "$APPS/helixer-docker_helixer_v0.3.0_cuda_11.2.0-cudnn8.sif" \
  Helixer.py \
  --species solanum_verrucosum \
  --fasta-path results/final_assembly/final_assembly.fa \
  --gff-output-path results/helixer/helixer.gff \
  --model-filepath "$APPS/helixer_model/land_plant_v0.3_a_0080.h5"

$AGAT_SP_EXTRACT_SEQUENCES -g results/helixer/helixer.gff -f results/final_assembly/final_assembly.fa -p -t cds -o results/helixer/helixer.pep.fa
$AGAT_SP_EXTRACT_SEQUENCES -g results/helixer/helixer.gff -f results/final_assembly/final_assembly.fa -t cds -o results/helixer/helixer.cds.fa