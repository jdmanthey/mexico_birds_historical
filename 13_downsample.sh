#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=downsample
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-208

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/downsample_helper.txt | cut -f1 | tail -n1 )
prop_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/downsample_helper.txt | cut -f3 | tail -n1 )


# downsample if prop_array < 1, copy if prop_array = 1       
if [[ $prop_array == "1" ]]; then 
	
	cp ${workdir}/01_bam_files/${basename_array}_final.bam ${workdir}/01b_bam_files/${basename_array}_final.bam

	cp ${workdir}/01_bam_files/${basename_array}_final.bam.bai ${workdir}/01b_bam_files/${basename_array}_final.bam.bai

else
	
	java -jar /home/jmanthey/picard.jar DownsampleSam \
	I=${workdir}/01_bam_files/${basename_array}_final.bam \
	O=${workdir}/01b_bam_files/${basename_array}_final.bam \
	P=${prop_array}
	
	# index the final bam file
	samtools index ${workdir}/01b_bam_files/${basename_array}_final.bam
	
fi






