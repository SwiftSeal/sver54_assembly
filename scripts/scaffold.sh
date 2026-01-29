#!/bin/bash

#SBATCH -p long
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -o logs/scaffold.%j.log
#SBATCH -e logs/scaffold.%j.log

set -euo pipefail

ASSEMBLY_PATH="results/hifiasm/hifiasm.p_ctg.fa"
READ_PATH="/mnt/shared/projects/jhi/potato/202210_Sver-HiFi_Moray/hic"
R1_1="$READ_PATH/sver_54_1_S1_R1_001.fastq.gz"
R1_2="$READ_PATH/sver_54_2_S2_R1_001.fastq.gz"
R2_1="$READ_PATH/sver_54_1_S1_R2_001.fastq.gz"
R2_2="$READ_PATH/sver_54_2_S2_R2_001.fastq.gz"

mkdir -p results/scaffolding

singularity exec -B /mnt/:/mnt/ containers/bwa.sif bwa index $ASSEMBLY_PATH

singularity exec -B /mnt/:/mnt/ containers/bwa.sif \
  bwa mem \
  -5SP \
  -T0 \
  -t $SLURM_CPUS_PER_TASK \
  $ASSEMBLY_PATH \
  <(zcat $R1_1 $R1_2) <(zcat $R2_1 $R2_2) \
  -o $TMPDIR/aligned.sam

singularity exec -B /mnt/:/mnt/ containers/pairtools.sif \
  pairtools parse \
  --min-mapq 40 \
  --walks-policy 5unique \
  --max-inter-align-gap 30 \
  --chroms-path $ASSEMBLY_PATH \
  --output results/scaffolding/parsed.pairsam \
  $TMPDIR/aligned.sam

singularity exec -B /mnt/:/mnt/ containers/pairtools.sif \
  pairtools sort \
  --nproc $SLURM_CPUS_PER_TASK \
  --output results/scaffolding/sorted.pairsam \
  results/scaffolding/parsed.pairsam

singularity exec -B /mnt/:/mnt/ containers/pairtools.sif \
  pairtools dedup \
  --mark-dups \
  --output-stats results/scaffolding/pairtools_stats.txt \
  --output results/scaffolding/dedup.pairsam \
  results/scaffolding/sorted.pairsam

singularity exec -B /mnt/:/mnt/ containers/pairtools.sif \
  pairtools split \
  --output-pairs results/scaffolding/mapped.pairs \
  --output-sam results/scaffolding/unsorted.pairsam \
  results/scaffolding/dedup.pairsam

singularity exec -B /mnt/:/mnt/ containers/samtools.sif \
  samtools sort \
  -@ $SLURM_CPUS_PER_TASK \
  -o results/scaffolding/mapped_pairtools.bam \
  results/scaffolding/unsorted.pairsam