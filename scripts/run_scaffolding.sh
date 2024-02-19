#!/bin/bash

#SBATCH -p long
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/scaffolding.%j.out
#SBATCH -e logs/scaffolding.%j.err

zcat "/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_1_S1_R1_001.fastq.gz" \
"/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_2_S2_R1_001.fastq.gz" \
> $TMPDIR/R1.fastq.gz

zcat "/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_1_S1_R2_001.fastq.gz" \
"/mnt/shared/projects/jhi/misc/miseq/NextSeq2000/230616_VH00754_44_AACTK5TM5/Analysis/1/Data/fastq/sver_54_2_S2_R2_001.fastq.gz" \
> $TMPDIR/R2.fastq.gz

source activate hic

bwa index \
assembly/hifiasm.p_ctg.fa

bwa mem \
-5SP \
-T0 \
-t 16 \
assembly/hifiasm.p_ctg.fa \
$TMPDIR/R1.fastq.gz $TMPDIR/R2.fastq.gz \
-o $TMPDIR/aligned.sam

rm $TMPDIR/R1.fastq.gz
rm $TMPDIR/R2.fastq.gz

pairtools parse \
--min-mapq 40 \
--walks-policy 5unique \
--max-inter-align-gap 30 \
--nproc-in 8 \
--nproc-out 8 \
--chroms-path assembly/hifiasm.p_ctg.fa \
$TMPDIR/aligned.sam > $TMPDIR/parsed.pairsam

pairtools sort \
--nproc 16 \
$TMPDIR/parsed.pairsam > $TMPDIR/sorted.pairsam

pairtools dedup \
--nproc-in 8 \
--nproc-out 8 \
--mark-dups \
--output-stats scaffolding/pairtools_stats.txt \
--output $TMPDIR/dedup.pairsam \
$TMPDIR/sorted.pairsam

pairtools split \
--nproc-in 8 \
--nproc-out 8 \
--output-pairs $TMPDIR/mapped.pairs \
--output-sam $TMPDIR/unsorted.sam \
$TMPDIR/dedup.pairsam

samtools sort -@ 16 -o scaffolding/mapped_pairtools.sam $TMPDIR/unsorted.sam

