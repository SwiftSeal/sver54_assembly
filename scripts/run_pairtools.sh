#!/bin/bash

#SBATCH -p long
#SBATCH -c 16
#SBATCH --exclude=n17-28-256-starbuck,n19-28-384-nicknack,n19-28-384-oddjob
#SBATCH --mem=48G
#SBATCH --export=ALL
#SBATCH -o logs/pairtools.%j.out
#SBATCH -e logs/pairtools.%j.err

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)
READ_PATH="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/hic"

source activate pairtools

bwa-mem2 index $ASSEMBLY_PATH

zcat "$READ_PATH/sver_54_1_S1_R1_001.fastq.gz" \
"$READ_PATH/sver_54_2_S2_R1_001.fastq.gz" \
> $TMPDIR/R1.fastq.gz

zcat "$READ_PATH/sver_54_1_S1_R2_001.fastq.gz" \
"$READ_PATH/sver_54_2_S2_R2_001.fastq.gz" \
> $TMPDIR/R2.fastq.gz

bwa-mem2 mem -5SP -T0 -t 16 $ASSEMBLY_PATH $TMPDIR/R1.fastq.gz $TMPDIR/R2.fastq.gz | samtools sort -o results/scaffolding/aligned.bam

pairtools parse \
  --min-mapq 40 \
  --walks-policy 5unique \
  --max-inter-align-gap 30 \
  --chroms-path $ASSEMBLY_PATH \
  --output results/scaffolding/parsed.pairsam.gz \
  results/scaffolding/aligned.bam

pairtools sort \
  --nproc 16 \
  --output results/scaffolding/sorted.pairsam.gz \
  results/scaffolding/parsed.pairsam.gz

pairtools dedup \
  --mark-dups \
  --output-stats results/scaffolding/pairtools_stats.txt \
  --output results/scaffolding/dedup.pairsam.gz \
  results/scaffolding/sorted.pairsam.gz

pairtools split \
  --output-pairs results/scaffolding/mapped.pairs \
  --output-sam results/scaffolding/unsorted.pairsam \
  results/scaffolding/dedup.pairsam.gz

samtools sort \
  -@ 16 \
  -o results/scaffolding/mapped_pairtools.bam \
  results/scaffolding/unsorted.pairsam
