#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/busco.%j.log
#SBATCH -e logs/busco.%j.log

mkdir $TMPDIR/fastas

cp results/hifiasm/hifiasm.p_ctg.fa $TMPDIR/fastas/
cp results/flye/assembly.fasta $TMPDIR/fastas/
cp results/assembly/cpc54.assembly.fa $TMPDIR/fastas/

singularity exec -B /mnt/:/mnt/ containers/busco.sif \
  busco \
  --force \
  -c $SLURM_CPUS_PER_TASK \
  -i $TMPDIR/fastas \
  --out_path results/busco/assembly \
  -m genome \
  -l solanales_odb12