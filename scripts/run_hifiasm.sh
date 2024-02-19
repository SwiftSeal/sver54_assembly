#!/bin/bash

#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=64gb
#SBATCH --export=ALL

HIFI_READS="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/Solanum_verrucosum_54_hifi_reads.fastq.gz"

source activate hifiasm

hifiasm --primary -t 32 $HIFI_READS -o assembly/hifiasm

awk '/^S/{{print ">"$2;print $3}}' assembly/hifiasm.p_ctg.gfa > assembly/hifiasm.p_ctg.fa

