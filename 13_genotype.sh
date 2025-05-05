#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=genotype
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=21
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-208

source activate bcftools

threads=20

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/genotype_helper.txt | tail -n1 | cut -f 1 )
reference=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/genotype_helper.txt | tail -n1 | cut -f 2 )

# define the reference genome
refgenome=/home/jmanthey/references/${reference}

# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads ${threads} \
-f ${refgenome} ${workdir}/01_bam_files/${basename_array}_final.bam | \
bcftools call -m --threads ${threads} -o ${workdir}/02_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/02_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/02_vcf/${basename_array}.vcf.gz

# filter individual vcf files
bcftools view -i 'MIN(DP)>7' ${workdir}/02_vcf/${basename_array}.vcf.gz > \
${workdir}/03_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/03_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/03_vcf/${basename_array}.vcf.gz

# genotype stats
echo ${basename_array} > ${basename_array}.stats

# number of genotyped sites passing minimum depth filter
echo "sites genotyped" >> ${basename_array}.stats
gzip -cd ${workdir}/03_vcf/${basename_array}.vcf.gz | grep -v "^#" | wc -l >> ${basename_array}.stats













