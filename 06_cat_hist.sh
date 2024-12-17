#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=cat
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-81

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_historical.txt | tail -n1 )

# cat all split files together for each individual
for i in $( ls ${workdir}/01_split_h/${basename_array}*dedup_R1* ); do
	cat $i >> ${workdir}/01_working_h/${basename_array}_R1.fastq.gz;
done

for i in $( ls ${workdir}/01_split_h/${basename_array}*dedup_R2* ); do
	cat $i >> ${workdir}/01_working_h/${basename_array}_R2.fastq.gz;
done
