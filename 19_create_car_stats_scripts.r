options(scipen=999)
project_directory <- "/lustre/scratch/jmanthey/08_mexico/06_stats_car"
directory_name <- "36_stats_car"
cluster <- "nocona"
max_number_jobs <- 400

# command to open environment you will need
env_command <- "source activate bcftools" # if none, just put nothing in the quotes

# are you calculating any stats (probably yes)
stat_question <- TRUE

# are you running a RAxML phylogeny?
phylo_question <- TRUE

# vcf header location or file to extract the header (should be gzipped)
vcf_header <- paste0(project_directory, "/CM027536.1.recode.vcf.gz")

# name of script to write
output_script_name <- "01_stat_array.sh"

# read in reference index
# filtered to only include genotyped chromosomes
ref_index <- read.table("GCA_001746935.2_mywa_2.1_genomic.fna.fai", stringsAsFactors=F)

# define window size
window_size <- 100000

######################################
######################################
## only edit above this (unless you know specific edits you need below)
######################################
######################################

# make directories
dir.create(directory_name)

# define intervals and write to helper files
stat_helper1 <- list()
stat_helper2 <- list()
stat_helper3 <- list()
counter <- 1
for(a in 1:nrow(ref_index)) {
  a_start <- 1
  a_end <- a_start + window_size - 1
  a_max <- ref_index[a,2]
  a_windows <- ceiling((a_max - a_start) / window_size)
  a_chromosome <- ref_index[a,1]
  
  # loop for defining helper info for each window
  for(b in 1:a_windows) {
    if(b == a_windows) {
      a_end <- a_max
    }
    stat_helper1[[counter]] <- a_chromosome
    stat_helper2[[counter]] <- a_start
    stat_helper3[[counter]] <- a_end
    
    a_start <- a_start + window_size
    a_end <- a_end + window_size
    counter <- counter + 1
  }
}
stat_helper1 <- unlist(stat_helper1)
stat_helper2 <- unlist(stat_helper2)
stat_helper3 <- unlist(stat_helper3)

# calculate number of array jobs
if(length(stat_helper3) > max_number_jobs) {
  n_jobs_per_array <- ceiling(length(stat_helper3) / max_number_jobs)
  n_array_jobs <- ceiling(length(stat_helper3) / n_jobs_per_array)
} else {
  n_array_jobs <- length(stat_helper3)
  n_jobs_per_array <- 1
}

stat_helper1 <- c(stat_helper1, rep("x", n_jobs_per_array - length(stat_helper3) %% n_jobs_per_array))
stat_helper2 <- c(stat_helper2, rep(1, n_jobs_per_array - length(stat_helper3) %% n_jobs_per_array))
stat_helper3 <- c(stat_helper3, rep(1, n_jobs_per_array - length(stat_helper3) %% n_jobs_per_array))
length(stat_helper3)
stat_helper <- data.frame(chrom=as.character(stat_helper1), start=as.numeric(stat_helper2), end=as.numeric(stat_helper3))
write.table(stat_helper, file=paste(directory_name, "/stat_helper.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)

# write the array script
a.script <- paste(directory_name, "/", output_script_name, sep="")
write("#!/bin/sh", file=a.script)
write("#SBATCH --chdir=./", file=a.script, append=T)
write(paste("#SBATCH --job-name=", "stats", sep=""), file=a.script, append=T)
write("#SBATCH --nodes=1 --ntasks=2", file=a.script, append=T)
write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
write("#SBATCH --time=48:00:00", file=a.script, append=T)
write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
write(paste("#SBATCH --array=1-", n_array_jobs, sep=""), file=a.script, append=T)
write("", file=a.script, append=T)
write(env_command, file=a.script, append=T)
write("", file=a.script, append=T)

write("# Set the number of runs that each SLURM task should do", file=a.script, append=T)
write(paste("PER_TASK=", n_jobs_per_array, sep=""), file=a.script, append=T)
write("", file=a.script, append=T)

write("# Calculate the starting and ending values for this task based", file=a.script, append=T)
write("# on the SLURM task and the number of runs per task.", file=a.script, append=T)
write("START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))", file=a.script, append=T)
write("END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))", file=a.script, append=T)
write("", file=a.script, append=T)

write("# Print the task and run range", file=a.script, append=T)
write("echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM", file=a.script, append=T)
write("", file=a.script, append=T)

write("# Run the loop of runs for this task.", file=a.script, append=T)	
write("for (( run=$START_NUM; run<=$END_NUM; run++ )); do", file=a.script, append=T)
write("\techo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run", file=a.script, append=T)
write("", file=a.script, append=T)

write("\tchrom_array=$( head -n${run} stat_helper.txt | cut -f1 | tail -n1 )", file=a.script, append=T)
write("", file=a.script, append=T)
write("\tstart_array=$( head -n${run} stat_helper.txt | cut -f2 | tail -n1 )", file=a.script, append=T)
write("", file=a.script, append=T)
write("\tend_array=$( head -n${run} stat_helper.txt | cut -f3 | tail -n1 )", file=a.script, append=T)
write("", file=a.script, append=T)

# add header to output file
header <- paste('\tgunzip -cd ', vcf_header, ' | grep "#" > ', project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf", sep="")
write(header, file=a.script, append=T)
write("", file=a.script, append=T)

#tabix command
tabix_command <- paste("\ttabix ", project_directory, "/${chrom_array}.recode.vcf.gz ${chrom_array}:${start_array}-${end_array} >> ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf", sep="")
write(tabix_command, file=a.script, append=T)
write("", file=a.script, append=T)




# Rscript command for stats
if(stat_question == T) {
	rscript_command <- paste("\tRscript _window_stat_calculations.r ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf popmap_stats.txt", sep="")
	write(rscript_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	# remove unnecessary files at end
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf", sep=""), file=a.script, append=T)
}

# raxml commands
if(phylo_question == T) {
	raxml_command <- paste("\traxmlHPC-PTHREADS-SSE3 -T 2 -f a -x 50 -m GTRCAT -p 253 -N 2 -s ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.fasta -n ${chrom_array}__${start_array}__${end_array}.tre -w ", project_directory, "/windows/", sep="")
	write(raxml_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	# remove unnecessary files at end
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.fasta", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/RAxML_bestTree.${chrom_array}__${start_array}__${end_array}.tre", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/RAxML_bipartitionsBranchLabels.${chrom_array}__${start_array}__${end_array}.tre", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/RAxML_info.${chrom_array}__${start_array}__${end_array}.tre", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/RAxML_bootstrap.${chrom_array}__${start_array}__${end_array}.tre", sep=""), file=a.script, append=T)
}



write("", file=a.script, append=T)

# finish
write("done", file=a.script, append=T)
