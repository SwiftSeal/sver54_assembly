#!/bin/bash

#SBATCH -p himem
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/quickmerge.%j.out
#SBATCH -e logs/quickmerge.%j.err

source activate quickmerge

mkdir -p quickmerge

nucmer -l 100 -prefix results/quickmerge/quickmerge \
  results/flye/assembly.fasta \
  assembly/hifiasm.p_ctg.fa

delta-filter -r -q -l 10000 results/quickmerge/quickmerge.delta > results/quickmerge/quickmerge.rq.delta

cd results/quickmerge # because quickmerge has stupid outputs
quickmerge \
  -hco 5 \
  -c 1.5 \
  -l 2000000 \
  -ml 10000 \
  -p quickmerge \
  -d quickmerge.rq.delta \
  -q ../assembly/hifiasm.p_ctg.fa \
  -r ../flye/assembly.fasta
