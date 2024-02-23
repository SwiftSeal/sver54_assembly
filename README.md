# Solanum verrucosum clone 54 assembly

This repository contains all the code used in the *S. verrucosum* genome assembly project.
Here's a brief overview of the repository:

* `results` - computational results are stored in this directory
* `scripts` - `SLURM` scripts used to execute analysis jobs are contained here.
* `analysis` - code used in the analysis of this genome is stored here.
* `envs` - `conda` configuration files for the environments.
* `logs` - for storing `SLURM` output files.

## Assembly

The initial genome assembly is produced via `run_hifiasm.sh` which uses `hifiasm` and HiFi reads.
This is the scaffolded using Hi-C data.
First, a pairtools pipeline is used - `run_pairtools.sh`.
This is followed by `run_yahs.sh` to scaffold the genome with `yahs`.
The scaffold is manually curated externally with the Juicer GUI.
