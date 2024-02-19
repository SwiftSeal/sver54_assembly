#!/bin/bash

#SBATCH -p medium
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/oatk.%j.out
#SBATCH -e logs/oatk.%j.err

HIFI_READS="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/Solanum_verrucosum_54_hifi_reads.fastq.gz"
MITO_DB_PATH="/mnt/shared/scratch/msmith/apps/oatk/OatkDB/v20230921/embryophyta_mito.fam"
CHLO_DB_PATH="/mnt/shared/scratch/msmith/apps/oatk/OatkDB/v20230921/embryophyta_pltd.fam"

source activate oatk

$APPS/oatk/oatk -c 100 -t 8 -m $MITO_DB_PATH -p $CHLO_DB_PATH -o organelle/oatk $HIFI_READS

