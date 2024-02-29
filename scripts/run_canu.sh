#!/bin/bash

#SBATCH -p long
#SBATCH --export=ALL

source activate canu

canu \
  -p canu \
  -d results/canu \
  genomeSize=700m \
  -nanopore results/guppy/reads.fastq
