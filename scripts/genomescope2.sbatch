#!/bin/env -S pixi exec -c conda-forge -c bioconda --spec kmc==3.2.4 --spec genomescope2==2.0.1 -- bash

#SBATCH -p medium
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -o logs/genomescope2.%j.log
#SBATCH -e logs/genomescope2.%j.log

HIFI_READS="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/Solanum_verrucosum_54_hifi_reads.fastq.gz"
K=21

kmc -k$K -t$SLURM_CPUS_PER_TASK -m64 -ci1 -cs10000 $HIFI_READS $TMPDIR/reads $TMPDIR

kmc_tools transform $TMPDIR/reads histogram $TMPDIR/reads.histo -cx10000

genomescope2 -i $TMPDIR/reads.histo -o results/genomescope2 -k $K
