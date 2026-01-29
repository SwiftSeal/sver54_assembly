#!/bin/env bash

#SBATCH --array=1
#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --mem=32gb
#SBATCH --gpus=1
#SBATCH --export=ALL
#SBATCH -o logs/deepsignal.%j.%a.log
#SBATCH -e logs/deepsignal.%j.%a.log

export HDF5_PLUGIN_PATH=$APPS/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

set -euo pipefail

ONT_DIRS=(
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/2023-04-11_ver54_p5"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/2023-04-13_ver54_p5"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Ver54_2_200422"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Ver54_200422"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Sver54_seq"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Sver_repeat_MS"
)

ONT_DIR=${ONT_DIRS[$SLURM_ARRAY_TASK_ID-1]}

mkdir -p results/deepsignal

pixi exec -c conda-forge -c bioconda --spec "setuptools<81" --spec ont-fast5-api multi_to_single_fast5 \
  -i $ONT_DIR \
  -s $TMPDIR/fast5_single \
  --recursive \
  -t 16

$APPS/ont-guppy/bin/guppy_basecaller -r -i $TMPDIR/fast5_single -s $TMPDIR/guppy -c dna_r9.4.1_450bps_sup_plant.cfg -x 'auto'

cat $TMPDIR/guppy/pass/*.fastq > $TMPDIR/reads.fastq
cat $TMPDIR/guppy/sequencing_summary.txt > $TMPDIR/sequencing_summary.txt
rm -r $TMPDIR/guppy

pixi exec -c conda-forge -c bioconda --spec ont-tombo=1.5.1 tombo preprocess annotate_raw_with_fastqs \
  --processes $SLURM_CPUS_PER_TASK \
  --overwrite \
  --fast5-basedir $TMPDIR/fast5_single \
  --fastq-filenames $TMPDIR/reads.fastq \
  --basecall-group Basecall_1D_000 \
  --basecall-subgroup BaseCalled_template \
  --sequencing-summary-filenames $TMPDIR/sequencing_summary.txt

pixi exec -c conda-forge -c bioconda --spec ont-tombo=1.5.1 tombo resquiggle \
  $TMPDIR/fast5_single \
  results/assembly/cpc.assembly.fa \
  --processes $SLURM_CPUS_PER_TASK \
  --overwrite \
  --corrected-group RawGenomeCorrected_000 \
  --basecall-group Basecall_1D_000

pixi exec -c conda-forge -c bioconda --spec deepsignal-plant=0.1.6 deepsignal_plant call_mods \
  --input_path $TMPDIR/fast5_single \
  --model_path $APPS/model.dp2.CNN.arabnrice2-1_120m_R9.4plus_tem.bn13_sn16.both_bilstm.epoch6.ckpt \
  --result_file "results/deepsignal/$(basename $ONT_DIR).tsv" \
  --corrected_group RawGenomeCorrected_000 \
  --motifs C \
  --nproc $SLURM_CPUS_PER_TASK \
  --nproc_gpu 4
