#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter_cat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-40

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds_thrush.txt | tail -n1 )

# filter for structure (minimum 53/56 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Cat_${region_array}_6.vcf --keep keep_cat.txt --thin 10000 \
--max-missing 0.94 --mac 4 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_pca/Cat_structure_${region_array}

# invariant and variant sites for stats from three different coverage datasets
# filter for structure (minimum 51/56 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Cat_${region_array}_6.vcf --keep keep_cat.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats6_cat/${region_array}

# invariant and variant sites for stats from three different coverage datasets
# filter for structure (minimum 51/56 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Cat_${region_array}_8.vcf --keep keep_cat.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats8_cat/${region_array}

# invariant and variant sites for stats from three different coverage datasets
# filter for structure (minimum 51/56 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Cat_${region_array}_10.vcf --keep keep_cat.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats10_cat/${region_array}


# bgzip and tabix index files that will be subdivided into windows
bgzip ${workdir}/06_stats6_cat/${region_array}.recode.vcf

bgzip ${workdir}/06_stats8_cat/${region_array}.recode.vcf

bgzip ${workdir}/06_stats10_cat/${region_array}.recode.vcf

#tabix
tabix -p vcf ${workdir}/06_stats6_cat/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/06_stats8_cat/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/06_stats10_cat/${region_array}.recode.vcf.gz

