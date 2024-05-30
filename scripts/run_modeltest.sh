#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/mafft:7.525--h031d066_1 mafft \
    --thread 8 \
    nbarc.fasta \
    > nbarc.aln.fa

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/modeltest-ng:0.1.7--h103dbdd_2 modeltest-ng \
  -i nbarc.aln.fa \
  -d aa \
  -p 8 \
  -o modeltest-ng

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/raxml-ng:1.2.2--h6d1f11b_0 raxml-ng \
  --threads 8 \
  --msa nbarc.aln.fa \
  --model JTT+G4