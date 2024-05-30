#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o logs/busco.%j.out
#SBATCH -e logs/busco.%j.err

BUSCO="singularity exec -B /mnt/:/mnt/ docker://ezlabgva/busco:v5.7.0_cv1 busco"

$BUSCO \
  -c 8 \
  -i results/braker/braker.pep.fa \
  -o braker \
  --out_path results/busco \
  -m prot \
  -l viridiplantae_odb10

$BUSCO \
  -c 8 \
  -i results/helixer/helixer.pep.fa \
  -o helixer \
  --out_path results/busco \
  -m prot \
  -l viridiplantae_odb10

$BUSCO \
  -c 8 \
  -i results/stringtie/transcripts.fa.transdecoder.pep \
  -o stringtie \
  --out_path results/busco \
  -m prot \
  -l viridiplantae_odb10

$BUSCO \
  -c 1 \
  -i results/final_annotation/final_annotation.pep.fa \
  -o final \
  --out_path results/busco \
  -m prot \
  -l viridiplantae_odb10
