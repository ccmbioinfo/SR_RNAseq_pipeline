#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --tmp=350G
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --constraint=AlmaLinux9

#module load python/3.12.4
module load Singularity
module load sqlite/3.47.0

source /hpf/largeprojects/tgnode/sandbox/rnaseq_repo/SR_RNAseq_pipeline_env/miniforge3/etc/profile.d/conda.sh
conda activate new_SR_RNAseq_pipeline

#source ${pipeline_location}/python_venv/bin/activate

snakemake --snakefile ${pipeline_location}/Snakefile --cores 32 \
--configfile $configfile \
--use-singularity \
--singularity-args "-B /hpf/largeprojects/,/hpf/projects/,$SCRATCH:$SCRATCH" \
--nolock --rerun-triggers mtime 
#--rerun-incomplete
#--rerun-triggers mtime
