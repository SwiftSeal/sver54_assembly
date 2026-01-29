#!/bin/env -S pixi exec -c conda-forge -c bioconda --spec samtools==1.20 -- bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH --export ALL
#SBATCH -o logs/sizes.%j.out
#SBATCH -e logs/sizes.%j.out

ASSEMBLY_FASTA=results/assembly/cpc54.assembly.fa

samtools faidx $ASSEMBLY_FASTA
cut -f1,2 $ASSEMBLY_FASTA.fai > $ASSEMBLY_FASTA.sizes

