#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=cat
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-142

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/base_rename.txt | tail -n1 )

cat ${workdir}/00_fastq_m/temp/${basename}_*R1_* >> ${workdir}/00_fastq_m/${basename}_R1.fastq.gz

cat ${workdir}/00_fastq_m/temp/${basename}_*R2_* >> ${workdir}/00_fastq_m/${basename}_R2.fastq.gz

