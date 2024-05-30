#!/bin/bash

#SBATCH -p himem
#SBATCH -c 1
#SBATCH --mem=200gb
#SBATCH --export=ALL
#SBATCH -o logs/deepsignal_freq.%j.out
#SBATCH -e logs/deepsignal_freq.%j.err

source activate deepsignal

cat results/deepsignal/* > $TMPDIR/mods.tsv

deepsignal_plant call_freq \
  --input_path $TMPDIR/mods.tsv \
  --result_file results/deepsignal/freq.tsv

python scripts/split_freq_file_by_5mC_motif.py --freqfile $TMPDIR/freq.tsv

# convert to bedgraph

awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$10}' results/deepsignal/freq.CG.tsv > $TMPDIR/freq.CG.bedgraph
sort -k1,1 -k2,2n $TMPDIR/freq.CG.bedgraph > results/deepsignal/freq.CG.bedgraph
awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$10}' results/deepsignal/freq.CHG.tsv > $TMPDIR/freq.CHG.bedgraph
sort -k1,1 -k2,2n $TMPDIR/freq.CHG.bedgraph > results/deepsignal/freq.CHG.bedgraph
awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$10}' results/deepsignal/freq.CHH.tsv > $TMPDIR/freq.CHH.bedgraph
sort -k1,1 -k2,2n $TMPDIR/freq.CHH.bedgraph > results/deepsignal/freq.CHH.bedgraph

ASSEMBLY_FASTA="results/final_assembly/final_assembly.fa"
$SAMTOOLS faidx "$ASSEMBLY_FASTA"
cut -f1,2 "$ASSEMBLY_FASTA.fai" > "$TMPDIR/chrom.sizes"

bedGraphToBigWig results/deepsignal/freq.CG.bedgraph $TMPDIR/chrom.sizes results/deepsignal/freq.CG.bw
bedGraphToBigWig results/deepsignal/freq.CHG.bedgraph $TMPDIR/chrom.sizes results/deepsignal/freq.CHG.bw
bedGraphToBigWig results/deepsignal/freq.CHH.bedgraph $TMPDIR/chrom.sizes results/deepsignal/freq.CHH.bw