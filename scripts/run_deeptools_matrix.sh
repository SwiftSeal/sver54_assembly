#!/bin/bash

#SBATCH -p himem
#SBATCH -c 1
#SBATCH --mem=200gb
#SBATCH --export=ALL
#SBATCH -o logs/deepsignal_freq.%j.out
#SBATCH -e logs/deepsignal_freq.%j.err

split -l 5000 results/final_annotation/final_annotation.longest.bed $TMPDIR/annotation_chunks

for chunk in $TMPDIR/annotation_chunks*; do
  echo "$chunk"
  singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 computeMatrix scale-regions \
    -S results/deepsignal/freq.CG.bw results/deepsignal/freq.CHG.bw results/deepsignal/freq.CHH.bw \
    -R $chunk \
    -m 2000 -a 2000 -b 2000 \
    -p 4 \
    -o ${chunk}.gz
done

cat $TMPDIR/*.gz | gunzip -c > results/deepsignal/gene_matrix.tab

awk '{print $1"\t"$4"\t"$5"\t"$9"\t"$7}' results/earlgrey/solanum_verrucosum_EarlGrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.filteredRepeats.gff > $TMPDIR/te.bed

split -l 5000 $TMPDIR/te.bed $TMPDIR/te_chunks

for chunk in $TMPDIR/te_chunks*; do
  echo "$chunk"
  singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 computeMatrix scale-regions \
    -S results/deepsignal/freq.CG.bw results/deepsignal/freq.CHG.bw results/deepsignal/freq.CHH.bw \
    -R $chunk \
    -m 2000 -a 2000 -b 2000 \
    -p 4 \
    -o ${chunk}.gz
done

cat $TMPDIR/*.gz | gunzip -c > results/deepsignal/te_matrix.tab