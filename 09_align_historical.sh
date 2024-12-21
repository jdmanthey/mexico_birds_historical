#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-81

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/align_helper.txt | tail -n1 | cut -f 1 )
reference=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/align_helper.txt | tail -n1 | cut -f 2 )

# define the reference genome
refgenome=/home/jmanthey/references/${reference}

# run bwa mem for unmerged reads
bwa-mem2 mem -t 12 ${refgenome} ${workdir}/01_cleaned_h/${basename_array}_R1.fastq.gz \
${workdir}/01_cleaned_h/${basename_array}_R2.fastq.gz > ${workdir}/01_bam_files/${basename_array}_unmerged.sam

# run bwa mem for merged reads
bwa-mem2 mem -t 12 ${refgenome} ${workdir}/01_cleaned_h/${basename_array}_U.fastq.gz \
> ${workdir}/01_bam_files/${basename_array}_merged.sam

# convert sam to bam for unmerged reads
samtools view -b -S -o ${workdir}/01_bam_files/${basename_array}_unmerged.bam \
${workdir}/01_bam_files/${basename_array}_unmerged.sam

# convert sam to bam for merged reads
samtools view -b -S -o \
${workdir}/01_bam_files/${basename_array}_merged.bam ${workdir}/01_bam_files/${basename_array}_merged.sam

# combine the two bam files
samtools cat ${workdir}/01_bam_files/${basename_array}_merged.bam \
${workdir}/01_bam_files/${basename_array}_unmerged.bam > ${workdir}/01_bam_files/${basename_array}.bam

# remove ungrouped sam and bam files
rm ${workdir}/01_bam_files/${basename_array}_unmerged.sam
rm ${workdir}/01_bam_files/${basename_array}_unmerged.bam
rm ${workdir}/01_bam_files/${basename_array}_merged.sam
rm ${workdir}/01_bam_files/${basename_array}_merged.bam

# clean up the bam file
java -jar /home/jmanthey/picard.jar CleanSam \
I=${workdir}/01_bam_files/${basename_array}.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
java -jar /home/jmanthey/picard.jar SortSam I=${workdir}/01_bam_files/${basename_array}_cleaned.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
java -jar /home/jmanthey/picard.jar AddOrReplaceReadGroups \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates prefinal bam file)
java -jar /home/jmanthey/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
O=${workdir}/01_bam_files/${basename_array}_prefinal.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the prefinal bam file
samtools index ${workdir}/01_bam_files/${basename_array}_prefinal.bam

source activate mapdamage

module load  gcc/10.1.0 gsl/2.7 r/4.3.0

# make mapdamage output directory
mkdir ${workdir}/01_bam_files/temp_${basename_array}

# run mapdamage to rescale bam based on damage estimates
mapDamage -i ${workdir}/01_bam_files/${basename_array}_prefinal.bam -r ${refgenome} \
-d ${workdir}/01_bam_files/temp_${basename_array} --merge-reference-sequences --rescale \
--rescale-out=${workdir}/01_bam_files/${basename_array}_final.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${basename_array}_final.bam

# remove the prefinal bam file and index
rm ${workdir}/01_bam_files/${basename_array}_prefinal.bam
rm ${workdir}/01_bam_files/${basename_array}_prefinal.bam.bai

