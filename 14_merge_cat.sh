#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge_cat
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-40

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds_thrush.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Cat_*vcf.gz > \
${workdir}/04_vcf/Cat_${region_array}.vcf

# filter for structure (minimum 53/56 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Cat_${region_array}.vcf \
--max-missing 0.94 --mac 2 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_pca/Cat_structure_${region_array}

