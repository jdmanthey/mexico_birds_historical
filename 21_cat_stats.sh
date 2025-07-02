#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=stat_cat
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-16

dir_array=$( head -n${SLURM_ARRAY_TASK_ID} stat_cat_helper.txt | tail -n1 | cut -f 1 )
header_array=$( head -n${SLURM_ARRAY_TASK_ID} stat_cat_helper.txt | tail -n1 | cut -f 2 )
out_array=$( head -n${SLURM_ARRAY_TASK_ID} stat_cat_helper.txt | tail -n1 | cut -f 3 )

# combine the output for different species into a single file each
grep 'pop1' ${dir_array}${header_array} > ${out_array}
for i in $( ls ${dir_array}*__stats.txt ); do grep -v 'pop1' $i >> ${out_array}; done

