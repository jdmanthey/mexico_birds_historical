#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge_myi
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-31

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds_warbler.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Myi_*_6.vcf.gz > \
${workdir}/04_vcf/Myi_${region_array}_6.vcf

bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Myi_*_8.vcf.gz > \
${workdir}/04_vcf/Myi_${region_array}_8.vcf

bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Myi_*_10.vcf.gz > \
${workdir}/04_vcf/Myi_${region_array}_10.vcf
