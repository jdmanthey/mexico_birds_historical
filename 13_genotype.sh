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

# full datasets
# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads ${threads} \
-f ${refgenome} ${workdir}/01_bam_files/${basename_array}_final.bam | \
bcftools call -m --threads ${threads} -o ${workdir}/02_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/02_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/02_vcf/${basename_array}.vcf.gz

# filter individual vcf files for different coverages
bcftools view -i 'MIN(DP)>5' ${workdir}/02_vcf/${basename_array}.vcf.gz > \
${workdir}/03_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/03_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/03_vcf/${basename_array}.vcf.gz


# downsampled datasets
# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads ${threads} \
-f ${refgenome} ${workdir}/01b_bam_files/${basename_array}_final.bam | \
bcftools call -m --threads ${threads} -o ${workdir}/02b_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/02b_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/02b_vcf/${basename_array}.vcf.gz

# filter individual vcf files for different coverages
bcftools view -i 'MIN(DP)>5' ${workdir}/02b_vcf/${basename_array}.vcf.gz > \
${workdir}/03b_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/03b_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/03b_vcf/${basename_array}.vcf.gz




