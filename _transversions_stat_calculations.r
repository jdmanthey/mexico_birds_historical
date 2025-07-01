
options(scipen=999)

# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# what do you want to calculate or use?
calc_heterozygosity <- TRUE		# calculate observed heterozygosity for every individual
calc_diff <- TRUE				# calculate pairwise FST and DXY for all pairwise comps. of pops with 
									# N >= calc_diff_minN


filename_simple <- args[1]
filename_simple <- strsplit(filename_simple, ".subset.vcf")[[1]][1]


# files in windows or whole chromosomes
vcf_type <- "windows" # two options, "windows" or "chromosomes"
chromosome_info <- c(0,0,0)


# determine base outname (edit as needed)
outname_base <- filename_simple


# minimum number of individuals for diff analyses 
calc_diff_minN <- 2

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# do not modify below this point
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################



####################################
####################################
####################################
####################################
# functions used in below sections
####################################
####################################
####################################
####################################


####################################
# read in a small vcf (don't use 
# for large vcf files)
####################################
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}


####################################
# observed heterozygosity
####################################
heterozygosity <- function(xxx, outname, chrom_info) {
	for(a in 10:ncol(xxx)) {
		# select this individual
		a_rep <- xxx[,a]
		# remove missing genotypes
    		a_rep <- a_rep[a_rep != "./."]
    		# count number of sites
    		a_total <- length(a_rep)
    		# find number of heterozygous sites
    		a_het <- length(a_rep[a_rep == "0/1"])
    	
    		# add to output
    		output_rep <- c(colnames(xxx[a]), "none", "heterozygosity", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), a_total, a_het, a_het/a_total)
    		write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
	}
}



####################################
# helper function to determine if a site 
# is only missing data for a 
# population (used across rows of vcf)
# used in differentiation function
####################################
# function to determine if a site is only missing data for a population (used across rows of vcf)
# used in differentiation function
total_missing <- function(xxx) {
	return(length(xxx[xxx != "./."]) > 0)
}	


####################################
# helper function to determine if a site 
# is polymorphic (to be applied 
# across rows) after removing missing
# used in differentiation function
####################################
polymorphic_function2 <- function(xxx) {
  xxx <- xxx[xxx != "./."]
  return(length(unique(xxx)))
}			


####################################
# calc differentiation stats
# fst is calculation of Reich et 
# al. 2009 for small sample sizes
# equation presented nicer in 
# Willing et al. 2012 page 9
# Dxy calculation is from Tavares 
# et al. 2018: pnas.1801832115
####################################
differentiation <- function(xxx1, xxx2, popname1, popname2, outname, chrom_info) {
	# remove sites that are completely missing from either population
	keep1 <- apply(xxx1[,4:ncol(xxx1)], 1, total_missing)
	keep2 <- apply(xxx2[,4:ncol(xxx2)], 1, total_missing)
	xxx1 <- xxx1[keep1 == TRUE & keep2 == TRUE, ]
	xxx2 <- xxx2[keep1 == TRUE & keep2 == TRUE, ]
  
	# count the total number of included genotyped sites at this point
	n_sites <- nrow(xxx1)
  
	# combine the two matrices to find sites that are variant w/in and between the two pops
	xxx_combined <- cbind(xxx1[,4:ncol(xxx1)], xxx2[,4:ncol(xxx2)])
	variant_sites <- as.numeric(apply(xxx_combined, 1, polymorphic_function2))
  
	# keep only variant sites
	xxx1_variant <- xxx1[variant_sites > 1, 4:ncol(xxx1)]
	xxx2_variant <- xxx2[variant_sites > 1, 4:ncol(xxx2)]
  
	# count the number of variant sites
	n_variant_sites <- nrow(xxx1_variant)
  
	# loop for each polymorphic site to calculate dxy
	dxy_all <- list()
	for(a in 1:nrow(xxx1_variant)) {
		a_rep1 <- as.character(xxx1_variant[a,])
		a_rep2 <- as.character(xxx2_variant[a,])
    
		# remove missing
		a_rep1 <- a_rep1[a_rep1 != "./."]
		a_rep2 <- a_rep2[a_rep2 != "./."]
    
		# measure proportion of reference allele 
		a_ref1 <- (length(a_rep1[a_rep1 == "0/0"]) * 2 + length(a_rep1[a_rep1 == "0/1"]) * 1) / (length(a_rep1) * 2)
		a_ref2 <- (length(a_rep2[a_rep2 == "0/0"]) * 2 + length(a_rep2[a_rep2 == "0/1"]) * 1) / (length(a_rep2) * 2)
    
		# calc dxy
		dxy_all[[a]] <- a_ref1 * (1 - a_ref2) + a_ref2 * (1 - a_ref1)
	}
	dxy_all <- sum(unlist(dxy_all)) / n_sites
  
  
	# loop for each polymorphic site to calculate fst
	numerator_all <- list()
	denominator_all <- list()
	for(a in 1:nrow(xxx1_variant)) {
		a_rep1 <- as.character(xxx1_variant[a,])
		a_rep2 <- as.character(xxx2_variant[a,])
    
		# remove missing
		a_rep1 <- a_rep1[a_rep1 != "./."]
		a_rep2 <- a_rep2[a_rep2 != "./."]
    
		# number of individuals per population
		pop1_ind_count <- length(a_rep1) 
		pop2_ind_count <- length(a_rep2)
    
		# non-reference allele counts
		alt_allele_count1 <- (2 * length(a_rep1[a_rep1 == "1/1"]) + 1 * length(a_rep1[a_rep1 == "0/1"]))
		alt_allele_count2 <- (2 * length(a_rep2[a_rep2 == "1/1"]) + 1 * length(a_rep2[a_rep2 == "0/1"]))
    
		# total allele counts
		all_allele_count1 <- 2 * length(a_rep1)
		all_allele_count2 <- 2 * length(a_rep2)
    
		# expected heterozygosity for each population
		expected_het1 <- (alt_allele_count1 * (all_allele_count1 - alt_allele_count1)) / 
							(all_allele_count1 * (all_allele_count1 - 1))
		expected_het2 <- (alt_allele_count2 * (all_allele_count2 - alt_allele_count2)) / 
							(all_allele_count2 * (all_allele_count2 - 1))
    
		# find the fst numerator and denominator values for this snp (they all get summed and divided for 
		# the final estimate)
		numerator_all[[a]] <- (alt_allele_count1 / (2 * pop1_ind_count) - 
								alt_allele_count2 / (2 * pop2_ind_count))^2 - 
								(expected_het1 / (2 * pop1_ind_count)) - 
      							(expected_het2 / (2 * pop2_ind_count))
		denominator_all[[a]] <- numerator_all[[a]] + expected_het1 + expected_het2		
	}
	# calculate total fst for this vcf
	fst_all <- sum(unlist(numerator_all)) / sum(unlist(denominator_all))
  
	# write to output for dxy and fst
	output_rep1 <- c(popname1, popname2, "Dxy", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_variant_sites, dxy_all)
    write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
    output_rep2 <- c(popname1, popname2, "Fst", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_variant_sites, fst_all)
    write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
}



####################################
####################################
####################################
####################################
# calculate stats, etc.
####################################
####################################
####################################
####################################

# read in input file
input_file <- read_vcf(args[1])

# read in populations
populations <- read.table(args[2], stringsAsFactors=F, header=T)

####################################
# make VCF file only have 
# genotypes in columns 10:ncol
####################################
for(a in 10:ncol(input_file)) {
	input_file[,a] <- substr(input_file[,a], 1, 3)
}


####################################
# remove sites that are not either 
# invariant or bi-allelic SNPs
####################################
input_file <- input_file[nchar(input_file$REF) == 1 & nchar(input_file$ALT) == 1, ]


####################################
# remove phasing information if 
# present
####################################
for(a in 10:ncol(input_file)) {
	input_file[,a] <- gsub("\\|", "/", input_file[,a])
}


####################################
# make stats output header if 
# any stats == TRUE
####################################
if(calc_heterozygosity == TRUE | calc_diff == TRUE ) {
	# define output name
	output_name <- paste(outname_base, "__stats.txt", sep="")
	# write output file
	write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=output_name, sep="\t")
}


####################################
# calc heterozygosity
####################################
if(calc_heterozygosity == TRUE) {
	# calculate heterozygosity for each individual
	heterozygosity(input_file, output_name, chromosome_info)
}



####################################
# calc Fst, Dxy
####################################
if(calc_diff == TRUE) {
	# modify popmap to only have populations with N >= calc_diff_minN
	pop_counts <- table(populations[,2])
	pop_counts <- names(pop_counts)[pop_counts >= calc_diff_minN]
	populations2 <- populations[populations[,2] %in% pop_counts,]
	
	all_combinations <- combn(unique(populations2[,2]), 2)
	
	for(a in 1:ncol(all_combinations)) {
		# define populations
		a_pop1 <- all_combinations[1,a]
		a_pop2 <- all_combinations[2,a]
	
		# subset vcf inputs
		a_input1 <- 	cbind(input_file[,c(1,4,5)], input_file[,colnames(input_file) 
			%in% populations2[populations2[,2] == a_pop1, 1]])
		a_input2 <- 	cbind(input_file[,c(1,4,5)], input_file[,colnames(input_file) 
			%in% populations2[populations2[,2] == a_pop2, 1]])

		# calculate stats
		differentiation(a_input1, a_input2, a_pop1, a_pop2, output_name, chromosome_info)
	}
}

