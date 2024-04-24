#!/bin/bash

#SBATCH -p long
#SBATCH -c 16
#SBATCH --exclude=n17-28-256-starbuck,n19-28-384-nicknack,n19-28-384-oddjob
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/pairtools.%j.out
#SBATCH -e logs/pairtools.%j.err

PAIRTOOLS="singularity exec -B /mnt/:/mnt/ https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0 pairtools"
SAMTOOLS="singularity exec -B /mnt/:/mnt/ https://depot.galaxyproject.org/singularity/samtools:1.19.1--h50ea8bc_0 samtools"

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)
READ_PATH="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/hic"

$PAIRTOOLS parse \
  --min-mapq 40 \
  --walks-policy 5unique \
  --max-inter-align-gap 30 \
  --chroms-path $ASSEMBLY_PATH \
  --output results/scaffolding/parsed.pairsam \
  results/scaffolding/aligned.sam

$PAIRTOOLS sort \
  --nproc 16 \
  --output results/scaffolding/sorted.pairsam \
  results/scaffolding/parsed.pairsam

$PAIRTOOLS dedup \
  --mark-dups \
  --output-stats results/scaffolding/pairtools_stats.txt \
  --output results/scaffolding/dedup.pairsam \
  results/scaffolding/sorted.pairsam

$PAIRTOOLS split \
  --output-pairs results/scaffolding/mapped.pairs \
  --output-sam results/scaffolding/unsorted.pairsam \
  results/scaffolding/dedup.pairsam

$SAMTOOLS sort \
  -@ 16 \
  -o results/scaffolding/mapped_pairtools.bam \
  results/scaffolding/unsorted.pairsam