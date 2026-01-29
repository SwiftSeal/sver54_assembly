#!/bin/env -S pixi exec hifiasm==0.24.0 -- bash

#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem 64gb
#SBATCH --export ALL

mkdir -p results/hifiasm/

HIFI_READS="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/Solanum_verrucosum_54_hifi_reads.fastq.gz"
K=21

hifiasm --primary -t $SLURM_CPUS_PER_TASK $HIFI_READS -o results/hifiasm/asm
