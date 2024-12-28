#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge_car
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
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Car_*vcf.gz > \
${workdir}/04_vcf/Car_${region_array}.vcf

# filter for structure (minimum 31/34 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Car_${region_array}.vcf \
--max-missing 0.9 --mac 2 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_pca/Car_structure_${region_array}

