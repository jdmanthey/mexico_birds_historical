#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=t_stats
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-4

source activate bcftools

vcf_array=$( head -n${SLURM_ARRAY_TASK_ID} ./t_stats_helper.txt | tail -n1 | cut -f 1 )
popmap_array=$( head -n${SLURM_ARRAY_TASK_ID} ./t_stats_helper.txt | tail -n1 | cut -f 2 )

Rscript _transversions_stat_calculations.r $vcf_array $popmap_array


