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