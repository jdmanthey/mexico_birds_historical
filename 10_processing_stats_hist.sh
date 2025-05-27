#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=stats
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-78

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_historical.txt | tail -n1 )

echo ${basename_array} > ${basename_array}.stats

# read1 preprocessing bp
/lustre/work/jmanthey/bbmap/stats.sh ${workdir}/00_fastq_h/${basename_array}_R1.fastq.gz | \
tail -n2 | head -n1 | cut -f5 | sed 's/,//g' | sed 's/\s//g' >> ${basename_array}.stats

# read2 preprocessing bp
/lustre/work/jmanthey/bbmap/stats.sh ${workdir}/00_fastq_h/${basename_array}_R2.fastq.gz | \
tail -n2 | head -n1 | cut -f5 | sed 's/,//g' | sed 's/\s//g' >> ${basename_array}.stats

# read 1 postprocessing bp
/lustre/work/jmanthey/bbmap/stats.sh ${workdir}/01_cleaned_h/${basename_array}_R1.fastq.gz | \
tail -n2 | head -n1 | cut -f5 | sed 's/,//g' | sed 's/\s//g' >> ${basename_array}.stats

# read 2 postprocessing bp
/lustre/work/jmanthey/bbmap/stats.sh ${workdir}/01_cleaned_h/${basename_array}_R2.fastq.gz | \
tail -n2 | head -n1 | cut -f5 | sed 's/,//g' | sed 's/\s//g' >> ${basename_array}.stats

# merged reads postprocessing bp
/lustre/work/jmanthey/bbmap/stats.sh ${workdir}/01_cleaned_h/${basename_array}_U.fastq.gz | \
tail -n2 | head -n1 | cut -f5 | sed 's/,//g' | sed 's/\s//g' >> ${basename_array}.stats

# samtools depth sum of aligned sites
samtools depth  ${workdir}/01_bam_files/${basename_array}_final.bam  |  \
awk '{sum+=$3} END { print sum }' >> ${basename_array}.stats
