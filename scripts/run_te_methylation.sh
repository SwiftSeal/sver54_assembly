#!/bin/bash

#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/te_methylation.%j.out
#SBATCH -e logs/te_methylation.%j.err

# Need to sort gff first
sort -k1,1 -k4,4n results/earlgrey/solanum_verrucosum_EarlGrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.filteredRepeats.gff > $TMPDIR/sorted.gff

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1 bedtools map \
    -a $TMPDIR/sorted.gff \
    -b results/deepsignal/freq.CG.bedgraph \
    -c 4 \
    -o mean,count \
    > results/deepsignal/earlgrey.CG.gff

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1 bedtools map \
    -a $TMPDIR/sorted.gff \
    -b results/deepsignal/freq.CHG.bedgraph \
    -c 4 \
    -o mean,count \
    > results/deepsignal/earlgrey.CHG.gff

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1 bedtools map \
    -a $TMPDIR/sorted.gff \
    -b results/deepsignal/freq.CHH.bedgraph \
    -c 4 \
    -o mean,count \
    > results/deepsignal/earlgrey.CHH.gff