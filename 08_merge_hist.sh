#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=seqprep
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-81

source activate hist_filter

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_historical.txt | tail -n1 )

# seqprep adapter trimming and merging short inserts that were sequenced
SeqPrep -f ${workdir}/01_working_h/${basename_array}_dedup_R1.fastq.gz \
-r ${workdir}/01_working_h/${basename_array}_dedup_R2.fastq.gz \
-1 ${workdir}/01_working_h/${basename_array}_unmerged_R1.fastq.gz \
-2 ${workdir}/01_working_h/${basename_array}_unmerged_R2.fastq.gz \
-A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
-y J \
-s ${workdir}/01_working_h/${basename_array}_merged.fastq.gz
