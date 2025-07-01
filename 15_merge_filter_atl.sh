#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge_atl
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-35

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds_sparrow.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/Atl_*.vcf.gz > \
${workdir}/04_vcf/Atl_${region_array}.vcf

bcftools merge -m id --regions ${region_array} ${workdir}/03b_vcf/Atl_*.vcf.gz > \
${workdir}/04b_vcf/Atl_${region_array}.vcf

# filter for structure 
vcftools --vcf ${workdir}/04_vcf/Atl_${region_array}.vcf --keep keep_atl.txt --thin 10000 \
--max-missing 0.8 --mac 4 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_pca/Atl_structure_${region_array}

vcftools --vcf ${workdir}/04b_vcf/Atl_${region_array}.vcf --keep keep_atl.txt --thin 10000 \
--max-missing 0.8 --mac 4 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05b_pca/Atl_structure_${region_array}

# invariant and variant sites for stats 
vcftools --vcf ${workdir}/04_vcf/Atl_${region_array}.vcf --keep keep_atl.txt \
--max-missing 0.65 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06_stats_atl/${region_array}

vcftools --vcf ${workdir}/04b_vcf/Atl_${region_array}.vcf --keep keep_atl.txt \
--max-missing 0.65 --max-alleles 2 \
--max-maf 0.49 --remove-indels --recode --recode-INFO-all \
--out ${workdir}/06b_stats_atl/${region_array}

# bgzip and tabix index files that will be subdivided into windows
bgzip ${workdir}/06_stats_atl/${region_array}.recode.vcf

bgzip ${workdir}/06b_stats_atl/${region_array}.recode.vcf

#tabix
tabix -p vcf ${workdir}/06_stats_atl/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/06b_stats_atl/${region_array}.recode.vcf.gz

# filter for biallelic sites
vcftools --vcf ${workdir}/04_vcf/Atl_${region_array}.vcf --keep keep_atl.txt \
--max-missing 0.65 --mac 2 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/07_transversions/Atl_${region_array}

vcftools --vcf ${workdir}/04b_vcf/Atl_${region_array}.vcf --keep keep_atl.txt \
--max-missing 0.65 --mac 2 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/07b_transversions/Atl_${region_array}

# extract header for transversions subset
grep "^#" ${workdir}/07_transversions/Atl_${region_array}.recode.vcf > \
${workdir}/07_transversions/Atl_${region_array}.subset.vcf

grep "^#" ${workdir}/07b_transversions/Atl_${region_array}.recode.vcf > \
${workdir}/07b_transversions/Atl_${region_array}.subset.vcf

# extract transversions
Rscript _extract_transversions.r ${workdir}/07_transversions/Atl_${region_array}.recode.vcf \
${workdir}/07_transversions/Atl_${region_array}.subset.vcf

Rscript _extract_transversions.r ${workdir}/07b_transversions/Atl_${region_array}.recode.vcf \
${workdir}/07b_transversions/Atl_${region_array}.subset.vcf

# remove intermediate unneeded vcfs 
rm ${workdir}/07_transversions/Atl_${region_array}.recode.vcf

rm ${workdir}/07b_transversions/Atl_${region_array}.recode.vcf







