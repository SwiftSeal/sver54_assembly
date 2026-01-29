#!/bin/bash

#SBATCH -p medium
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH -o logs/get_containers.%j.log
#SBATCH -e logs/get_containers.%j.log

mkdir -p containers

singularity pull containers/flye.sif docker://quay.io/biocontainers/flye:2.9.5--py312h5e9d817_2
singularity pull containers/hifiasm.sif https://depot.galaxyproject.org/singularity/hifiasm:0.24.0--h5ca1c30_0
singularity pull containers/kmc.sif docker://quay.io/biocontainers/kmc:3.2.4--h5ca1c30_4
singularity pull containers/genomescope2.sif docker://quay.io/biocontainers/genomescope2:2.0.1--py313r44hdfd78af_1
singularity pull containers/oatk.sif docker://quay.io/biocontainers/oatk:1.0--h577a1d6_1
singularity pull containers/bwa.sif https://depot.galaxyproject.org/singularity/bwa:0.7.16--pl5.22.0_0
singularity pull containers/pairtools.sif https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0
singularity pull containers/samtools.sif https://depot.galaxyproject.org/singularity/samtools:1.19.1--h50ea8bc_0
singularity pull containers/busco.sif docker://quay.io/biocontainers/busco:6.0.0--pyhdfd78af_0
singularity pull containers/edta.sif docker://quay.io/biocontainers/edta:2.2.2--hdfd78af_1
singularity pull containers/agat.sif docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0
singularity pull containers/helixer.sif docker://gglyptodon/helixer-docker:helixer_v0.3.0_cuda_11.2.0-cudnn8
singularity pull containers/minimap2.sif docker://quay.io/biocontainers/minimap2:2.30--h577a1d6_0
singularity pull containers/winnowmap.sif docker://quay.io/biocontainers/winnowmap:2.03--h5ca1c30_4
singularity pull containers/gci.sif docker://quay.io/biocontainers/gci:1.0--hdfd78af_0