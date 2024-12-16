#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=dedup
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=24G
#SBATCH --array=1-203

source activate hist_filter

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

read1=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/h_split.txt | tail -n1 | cut -f 1)
read2=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/h_split.txt | tail -n1 | cut -f 2)

read1=${read1%.fastq.gz}
read2=${read2%.fastq.gz}

# remove adaptase tails from IDT libraries
# also remove first couple bp from read 1 since there are many Ns
seqtk trimfq -b 10 ${workdir}/01_split_h/${read2}.fastq.gz | gzip > ${workdir}/01_split_h/${read2}_trimmed.fastq.gz

seqtk trimfq -b 2 ${workdir}/01_split_h/${read1}.fastq.gz | gzip > ${workdir}/01_split_h/${read1}_trimmed.fastq.gz

# remove duplicates
hts_SuperDeduper -1 ${workdir}/01_split_h/${read1}_trimmed.fastq.gz \
-2 ${workdir}/01_split_h/${read2}_trimmed.fastq.gz \
-L ${read2}_dedup.stats -f ${workdir}/01_split_h/${read2}_dedup
