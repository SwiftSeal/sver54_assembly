#!/bin/env -S pixi exec -c conda-forge -c bioconda --spec oatk==1.0  -- bash

#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH --export ALL
#SBATCH -o logs/oatk.%j.log
#SBATCH -e logs/oatk.%j.log

HIFI_READS="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/Solanum_verrucosum_54_hifi_reads.fastq.gz"
MITO_DB_PATH="results/organelle/OatkDB-1.0/v20230921/embryophyta_mito.fam"
CHLO_DB_PATH="results/organelle/OatkDB-1.0/v20230921/embryophyta_pltd.fam"

mkdir -p results/organelle
cd results/organelle
curl -LO https://github.com/c-zhou/OatkDB/archive/refs/tags/1.0.tar.gz
tar -xzf 1.0.tar.gz && rm 1.0.tar.gz
cd ../..


oatk -c 100 -t $SLURM_CPUS_PER_TASK -m $MITO_DB_PATH -p $CHLO_DB_PATH -o results/organelle/oatk $HIFI_READS

