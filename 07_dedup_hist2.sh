#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=dedup2
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=24G
#SBATCH --array=1-78

source activate hist_filter

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_historical.txt | tail -n1 )

# remove duplicates
hts_SuperDeduper -1 ${workdir}/01_working_h/${basename_array}_R1.fastq.gz \
-2 ${workdir}/01_working_h/${basename_array}_R2.fastq.gz \
-L ${basename_array}_dedup.stats -f ${workdir}/01_working_h/${basename_array}_dedup
