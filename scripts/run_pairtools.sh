#!/bin/bash

#SBATCH -p long
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/scaffolding.%j.out
#SBATCH -e logs/scaffolding.%j.err

ASSEMBLY_PATH=$(realpath results/quickmerge/merged_quickmerge.fasta)

source activate pairtools

cd $TMPDIR

bwa index $ASSEMBLY_PATH

zcat "/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_1_S1_R1_001.fastq.gz" \
"/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_2_S2_R1_001.fastq.gz" \
> R1.fastq.gz

zcat "/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_1_S1_R2_001.fastq.gz" \
"/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_2_S2_R2_001.fastq.gz" \
> R2.fastq.gz

bwa mem -5SP -T0 -t 16 $ASSEMBLY_PATH R1.fastq.gz R2.fastq.gz | samtools sort -o aligned.bam

rm R1.fastq.gz R2.fastq.gz

pairtools parse \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --chroms-path $ASSEMBLY_PATH \
    --output parsed.pairsam.gz \
    aligned.bam

rm aligned.bam

pairtools sort \
    --nproc 16 \
    --output sorted.pairsam.gz \
    parsed.pairsam.gz

rm parsed.pairsam.gz

pairtools dedup \
    --mark-dups \
    --output-stats pairtools_stats.txt \
    --output dedup.pairsam.gz \
    sorted.pairsam.gz

rm sorted.pairsam.gz

pairtools split \
    --output-pairs mapped.pairs \
    --output-sam unsorted.pairsam \
    dedup.pairsam

rm dedup.pairsam

samtools sort -@ 16 -o mapped_pairtools.bam $TMPDIR/unsorted.parsam

cd -

mkdir -p results/scaffolding
cp $TMPDIR/mapped_pairtools.bam results/scaffolding/
cp $TMPDIR/pairtools_stats.txt results/scaffolding/
