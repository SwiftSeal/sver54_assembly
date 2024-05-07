#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --mem=128gb
#SBATCH --gpus=2
#SBATCH --export=ALL
#SBATCH -o logs/deepsignal.%j.out
#SBATCH -e logs/deepsignal.%j.err

export HDF5_PLUGIN_PATH=/mnt/shared/scratch/msmith/apps/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

ONT_DIRS=(
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/2023-04-11_ver54_p5"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/2023-04-13_ver54_p5"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Ver54_2_200422"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Ver54_200422"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Sver54_seq"
"/mnt/shared/projects/jhi/potato/202205_Sver-ONT_Moray/Sver_repeat_MS"
)

mkdir -p results/deepsignal
cp results/final_assembly/final_assembly.fa $TMPDIR/genome.fa

source activate tombo

for ONT_DIR in ${ONT_DIRS[@]}; do
  multi_to_single_fast5 \
    -i $ONT_DIR \
    -s $TMPDIR/fast5_single \
    --recursive \
    -t 16
  
  $APPS/ont-guppy/bin/guppy_basecaller -r -i $TMPDIR/fast5_single -s $TMPDIR/guppy -c dna_r9.4.1_450bps_sup_plant.cfg -x 'auto'
  
  cat $TMPDIR/guppy/pass/*.fastq > $TMPDIR/reads.fastq
  cat $TMPDIR/guppy/sequencing_summary.txt > $TMPDIR/sequencing_summary.txt
  rm -r $TMPDIR/guppy
  
  tombo preprocess annotate_raw_with_fastqs \
    --processes 16 \
    --overwrite \
    --fast5-basedir $TMPDIR/fast5_single \
    --fastq-filenames $TMPDIR/reads.fastq \
    --basecall-group Basecall_1D_000 \
    --basecall-subgroup BaseCalled_template \
    --sequencing-summary-filenames $TMPDIR/sequencing_summary.txt
  
  tombo resquiggle \
    $TMPDIR/fast5_single \
    $TMPDIR/genome.fa \
    --processes 16 \
    --overwrite \
    --corrected-group RawGenomeCorrected_000 \
    --basecall-group Basecall_1D_000
  
  source deactivate
  source activate deepsignal
  
  deepsignal_plant call_mods --input_path $TMPDIR/fast5_single \
    --model_path $APPS/model.dp2.CNN.arabnrice2-1_120m_R9.4plus_tem.bn13_sn16.both_bilstm.epoch6.ckpt \
    --result_file "results/deepsignal/$(basename $ONT_DIR).tsv" \
    --corrected_group RawGenomeCorrected_000 \
    --motifs C \
    --nproc 16 \
    --nproc_gpu 4
done

cat results/deepsignal/*.tsv > results/deepsignal/mods.tsv

deepsignal_plant call_freq --input results/deepsignal/mods.tsv \
  --result_file results/deepsignal/freq.tsv