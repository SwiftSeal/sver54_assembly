#!/bin/bash

#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/merfin.%j.out
#SBATCH -e logs/merfin.%j.err

source activate merfin

ASSEMBLY_FASTA="/mnt/shared/scratch/msmith/solanum_verrucosum/results/genome/solanum_verrucosum.fa"
HIFI_READS="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/Solanum_verrucosum_54_hifi_reads.fastq.gz"

# count 21mers
meryl count k=21 "$HIFI_READS" output "$TMPDIR/hifi.meryl"
meryl count k=21 "$ASSEMBLY_FASTA" output "$TMPDIR/assembly.meryl"

# exclude kmers with count = 1 in the reads
meryl greater-than 1 "$TMPDIR/hifi.meryl" output "$TMPDIR/hifi.gt1.meryl"

merfin -completeness \
    -threads 16 \
    -seqmers "$TMPDIR/assembly.meryl" \
    -readmers "$TMPDIR/hifi.gt1.meryl" \
    -peak 48.8

merfin -hist \
    -threads 16 \
    -seqmers "$TMPDIR/assembly.meryl" \
    -readmers "$TMPDIR/hifi.gt1.meryl" \
    -peak 48.8 \
    -output "results/merfin.hist"
