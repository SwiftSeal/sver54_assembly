#!/bin/env -S pixi exec -c conda-forge -c bioconda --spec bowtie2==2.5.3 --spec samtools==1.20 --spec deeptools==3.5.5 -- bash

#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH --export=ALL
#SBATCH -o logs/cenh3.%j.out
#SBATCH -e logs/cenh3.%j.out

ACCESSION="SRR18548893"
OUTDIR="results/cenh3"

mkdir -p $OUTDIR

wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR185/093/$ACCESSION/$ACCESSION.fastq.gz

bowtie2-build results/assembly/cpc54.assembly.fa results/assembly/cpc54.assembly

bowtie2 -p $SLURM_CPUS_PER_TASK -x results/assembly/cpc54.assembly -U $OUTDIR/$ACCESSION.fastq.gz | \
  samtools sort -o $OUTDIR/$ACCESSION.bam -

samtools index $OUTDIR/$ACCESSION.bam

bamCoverage -b $OUTDIR/$ACCESSION.bam -o $OUTDIR/$ACCESSION.bw
