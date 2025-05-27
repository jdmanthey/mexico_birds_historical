#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=split
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-78

source activate hist_filter

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_historical.txt | tail -n1 )

seqkit split2 -1 ${workdir}/00_fastq_h/${basename_array}_R1.fastq.gz \
-2 ${workdir}/00_fastq_h/${basename_array}_R2.fastq.gz \
-s 100000000 -O ${workdir}/01_split_h/ -e .gz 
