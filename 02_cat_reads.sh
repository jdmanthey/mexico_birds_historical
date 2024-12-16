#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=cat
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-81

source activate hist_filter

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_historical.txt | tail -n1 )


cat ${workdir}/00_lane1_data/${basename_array}_R1.fastq.gz >> \
${workdir}/00_fastq_h/${basename_array}_R1.fastq.gz

cat ${workdir}/00_lane2_data/${basename_array}_R1.fastq.gz >> \
${workdir}/00_fastq_h/${basename_array}_R1.fastq.gz

cat ${workdir}/00_lane1_data/${basename_array}_R2.fastq.gz >> \
${workdir}/00_fastq_h/${basename_array}_R2.fastq.gz

cat ${workdir}/00_lane2_data/${basename_array}_R2.fastq.gz >> \
${workdir}/00_fastq_h/${basename_array}_R2.fastq.gz


