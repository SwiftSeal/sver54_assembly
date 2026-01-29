#!/bin/bash

#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH --export=ALL
#SBATCH -o logs/earlgrey.%j.log
#SBATCH -e logs/earlgrey.%j.log

source activate earlgrey

cp results/final_assembly/final_assembly.fa $TMPDIR
cd $TMPDIR

earlGrey -g final_assembly.fa -s solanum_verrucosum -o earlgrey -t 32

cp -r earlgrey $SLURM_SUBMIT_DIR/results/