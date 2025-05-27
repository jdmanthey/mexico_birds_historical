#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=align
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-130

source activate bcftools

threads=12

# define main working directory
workdir=/lustre/scratch/jmanthey/08_mexico

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/align_helper_modern.txt | tail -n1 | cut -f 1 )
reference=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/align_helper_modern.txt | tail -n1 | cut -f 2 )

# define the reference genome
refgenome=/home/jmanthey/references/${reference}

# run bbduk
/lustre/work/jmanthey/bbmap/bbduk.sh \
in1=${workdir}/00_fastq_m/${basename_array}_R1.fastq.gz \
in2=${workdir}/00_fastq_m/${basename_array}_R2.fastq.gz \
out1=${workdir}/01_cleaned_m/${basename_array}_R1.fastq.gz \
out2=${workdir}/01_cleaned_m/${basename_array}_R2.fastq.gz \
minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 \
ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bwa mem
bwa-mem2 mem -t ${threads} ${refgenome} \
${workdir}/01_cleaned_m/${basename_array}_R1.fastq.gz \
${workdir}/01_cleaned_m/${basename_array}_R2.fastq.gz > \
${workdir}/01_bam_files/${basename_array}.sam

# filter for mapped reads
samtools view -b -f 2 -@ ${threads} \
-o ${workdir}/01_bam_files/${basename_array}.bam \
${workdir}/01_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01_bam_files/${basename_array}.sam

# clean up the bam file
java -jar /home/jmanthey/picard.jar CleanSam \
I=${workdir}/01_bam_files/${basename_array}.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
java -jar /home/jmanthey/picard.jar SortSam \
I=${workdir}/01_bam_files/${basename_array}_cleaned.bam \
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

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
java -jar /home/jmanthey/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${basename_array}_final.bam










