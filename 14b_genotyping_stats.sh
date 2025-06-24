#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=gen_stat
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-208

source activate bcftools

threads=2

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/genotype_helper.txt | tail -n1 | cut -f 1 )

# alignment stats
echo ${basename_array} > ${basename_array}.stats

# samtools depth sum of aligned sites
echo "sum depth full" >> ${basename_array}.stats
samtools depth  ${workdir}/01_bam_files/${basename_array}_final.bam  |  awk '{sum+=$3} END { print sum }' >> ${basename_array}.stats

# samtools depth sum of aligned sites
echo "sum depth downsampled" >> ${basename_array}.stats
samtools depth  ${workdir}/01b_bam_files/${basename_array}_final.bam  |  awk '{sum+=$3} END { print sum }' >> ${basename_array}.stats

# number of genotyped sites passing minimum depth filter
echo "sites genotyped full" >> ${basename_array}.stats
gzip -cd ${workdir}/03_vcf/${basename_array}.vcf.gz | grep -v "^#" | wc -l >> ${basename_array}.stats

# number of genotyped sites passing minimum depth filter
echo "sites genotyped downsampled" >> ${basename_array}.stats
gzip -cd ${workdir}/03b_vcf/${basename_array}.vcf.gz | grep -v "^#" | wc -l >> ${basename_array}.stats



