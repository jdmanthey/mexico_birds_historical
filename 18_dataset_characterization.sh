#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=data_char
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-8

workdir=/lustre/scratch/jmanthey/08_mexico

dir_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/39_data_char/dir_helper.txt | tail -n1 )

cd ${workdir}/${dir_array} 

echo ${dir_array} 

for i in $( ls *gz ); do
	vcftools --gzvcf $i --mac 1 --out ${i%.recode.vcf.gz}
done
