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
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Car_*.vcf.gz > \
${workdir}/04_vcf/Car_${region_array}.vcf

bcftools merge -m id --regions ${region_array} ${workdir}/03b_vcf/Car_*.vcf.gz > \
${workdir}/04b_vcf/Car_${region_array}.vcf

# filter for structure (minimum 31/33 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Car_${region_array}.vcf --keep keep_car.txt --thin 10000 \
--max-missing 0.92 --mac 4 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_pca/Car_structure_${region_array}

vcftools --vcf ${workdir}/04b_vcf/Car_${region_array}.vcf --keep keep_car.txt --thin 10000 \
--max-missing 0.92 --mac 4 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05b_pca/Car_structure_${region_array}

# invariant and variant sites for stats 
vcftools --vcf ${workdir}/04_vcf/Car_${region_array}.vcf --keep keep_car.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats_car/${region_array}

vcftools --vcf ${workdir}/04b_vcf/Car_${region_array}.vcf --keep keep_car.txt \
--max-missing 0.9 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06b_stats_car/${region_array}

# bgzip and tabix index files that will be subdivided into windows
bgzip ${workdir}/06_stats_car/${region_array}.recode.vcf

bgzip ${workdir}/06b_stats_car/${region_array}.recode.vcf

#tabix
tabix -p vcf ${workdir}/06_stats_car/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/06b_stats_car/${region_array}.recode.vcf.gz

# filter for biallelic sites (minimum 31/33 individuals to keep a site)
vcftools --vcf ${workdir}/04_vcf/Car_${region_array}.vcf --keep keep_car.txt \
--max-missing 0.92 --mac 2 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/07_transversions/Car_${region_array}

vcftools --vcf ${workdir}/04b_vcf/Car_${region_array}.vcf --keep keep_car.txt \
--max-missing 0.92 --mac 2 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/07b_transversions/Car_${region_array}

# extract header for transversions subset
grep "^#" ${workdir}/07_transversions/Car_${region_array}.recode.vcf > \
${workdir}/07_transversions/Car_${region_array}.subset.vcf

grep "^#" ${workdir}/07b_transversions/Car_${region_array}.recode.vcf > \
${workdir}/07b_transversions/Car_${region_array}.subset.vcf

# extract transversions
Rscript _extract_transversions.r ${workdir}/07_transversions/Car_${region_array}.recode.vcf \
${workdir}/07_transversions/Car_${region_array}.subset.vcf

Rscript _extract_transversions.r ${workdir}/07b_transversions/Car_${region_array}.recode.vcf \
${workdir}/07b_transversions/Car_${region_array}.subset.vcf

# remove intermediate unneeded vcfs 
rm ${workdir}/07_transversions/Car_${region_array}.recode.vcf

rm ${workdir}/07b_transversions/Car_${region_array}.recode.vcf







