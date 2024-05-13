#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --gpus=1
#SBATCH --export=ALL

singularity exec --nv -B /mnt/:/mnt/ "$APPS/helixer-docker_helixer_v0.3.0_cuda_11.2.0-cudnn8.sif" \
  Helixer.py \
  --species solanum_verrucosum \
  --fasta-path results/final_assembly/final_assembly.fa \
  --gff-output-path helixer.gff \
  --model-filepath "$APPS/helixer_model/land_plant_v0.3_a_0080.h5"