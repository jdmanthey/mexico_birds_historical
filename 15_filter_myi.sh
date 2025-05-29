#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter_myi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-31

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds_warbler.txt | tail -n1 )

# filter for structure (minimum 42/44 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Myi_${region_array}_6.vcf --keep keep_myi.txt --thin 10000 \
--max-missing 0.95 --mac 4 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_pca/Myi_structure_${region_array}

# invariant and variant sites for stats from three different coverage datasets
# filter for structure (minimum 40/44 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Myi_${region_array}_6.vcf --keep keep_myi.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats6_myi/${region_array}

# invariant and variant sites for stats from three different coverage datasets
# filter for structure (minimum 40/44 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Myi_${region_array}_8.vcf --keep keep_myi.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats8_myi/${region_array}

# invariant and variant sites for stats from three different coverage datasets
# filter for structure (minimum 40/44 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Myi_${region_array}_10.vcf --keep keep_myi.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats10_myi/${region_array}


# bgzip and tabix index files that will be subdivided into windows
bgzip ${workdir}/06_stats6_myi/${region_array}.recode.vcf

bgzip ${workdir}/06_stats8_myi/${region_array}.recode.vcf

bgzip ${workdir}/06_stats10_myi/${region_array}.recode.vcf

#tabix
tabix -p vcf ${workdir}/06_stats6_myi/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/06_stats8_myi/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/06_stats10_myi/${region_array}.recode.vcf.gz

