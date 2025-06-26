options(scipen=999)

# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# what do you want to calculate or use?
calc_heterozygosity <- TRUE		# calculate observed heterozygosity for every individual
calc_pi <- TRUE					# calculate pi, theta, Tajima's D for each population with N >= calc_pi_minN
calc_diff <- TRUE				# calculate pairwise FST and DXY for all pairwise comps. of pops with 
									# N >= calc_diff_minN
calc_FIS <- FALSE				# calculate FIS for all pops with N >= calc_FIS_minN
calc_genetic_distance <- FALSE	# calculate genetic distance matrices
calc_polymorphisms <- FALSE		# calculate private, shared, fixed polymorphisms in all pops with 
									# N >= calc_polymorphisms_minN
calc_titv <- TRUE				# calculate transition / transversion ratio of entire dataset
make_FASTA <- TRUE				# make a FASTA file (e.g., for RAxML)


# files in windows or whole chromosomes
vcf_type <- "windows" # two options, "windows" or "chromosomes"

# base file name for parsing chromosome and window bounds if applicable
# modify this section as needed to work with your file names and filtering strategy
# the chromosome_info object should end as a vector with 3 strings -> chromosome, start, end
# if any of this info is missing, feel free to just make those values 0 so they are still a number, etc.
# modify the below lines as needed 
filename_simple <- args[1]
filename_simple <- strsplit(filename_simple, "windows/")[[1]][2]
filename_simple <- strsplit(filename_simple, ".recode")[[1]][1]
if(vcf_type == "chromosomes") {
	chromosome_info <- c(strsplit(filename_simple, ".recode.vcf")[[1]][1], 0, 0)
} else if( vcf_type == "windows") {
	chromosome_info <- c(strsplit(filename_simple, "__")[[1]][1], 
                    as.numeric(strsplit(strsplit(filename_simple, "__")[[1]][2], "__")[[1]][1]),
                    as.numeric(strsplit(filename_simple, "__")[[1]][3]))
}

# determine base outname (edit as needed)
outname_base <- strsplit(filename_simple, ".recode.vcf")[[1]][1]


# define minimum number of sites to keep a fasta file 
min_sites <- 10000


# minimum number of individuals for certain analyses (only used if those stats are calculated)
# recommended 3 or higher for each of the below, will not work with values of 0 or 1
# recommended 5 or higher for FIS
calc_pi_minN <- 2
calc_diff_minN <- 2
calc_polymorphisms_minN <- 3
calc_FIS_minN <- 5

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
# helper function to determine
# if a site is polymorphic to be
# applied across rows (e.g., a 
# subset of a VCF file)
####################################
# used in pi_tajima_theta function and fis function
polymorphic_function <- function(xxx) {
	xxx <- xxx[xxx != "./."]
	return(length(unique(xxx)))
}

####################################
# Fis helper function to get 
# expected het for each SNP to be
# applied across rows from VCF
# unbiased He from Hahn book 
# eqn. 3.1 
####################################
expected_heterozygosity <- function(xxx) {
	# remove missing genotypes
	xxx <- xxx[xxx != "./."]
    # sample size
	sample_size <- length(xxx) * 2
	# allele freq
	p <- (length(xxx[xxx == "0/0"]) * 2 + length(xxx[xxx == "0/1"])) / sample_size
	# expected heterozygosity = 2pq (with sample size correction)
	eh <- sample_size / (sample_size - 1) * (1 - p^2 - (p-1)^2)
	# output
	return(eh)
}


####################################
# calculate Fis
# Nei 1977 eqn. 14: F-statistics and 
# analysis of gene diversity in 
# subdivided populations
# He from previous function and 
# citation
####################################
fis <- function(xxx, popmap, outname, chrom_info) {
	for(a in 1:nrow(popmap)) {
		# select this individual
		a_rep <- xxx[,colnames(xxx) %in% popmap[a,1]]
		# select all individuals in the same population
		a_pop <- xxx[,colnames(xxx) %in% popmap[popmap[,2] == popmap[a,2],1]]
		# only keep polymorphic sites
		a_rep <- a_rep[as.numeric(apply(a_pop, 1, polymorphic_function)) > 1]
		a_pop <- a_pop[as.numeric(apply(a_pop, 1, polymorphic_function)) > 1,]		
		# calculate expected heterozygosity for each site
		a_eh <- apply(a_pop, 1, expected_heterozygosity)
		# average expected heterozygosity across all sites
		a_eh <- mean(na.omit(a_eh))
		# calculate observed heterozygosity for this individual
		# remove missing genotypes
		a_rep <- a_rep[a_rep != "./."]
		# find observed heterozygosity
		a_het <- length(a_rep[a_rep == "0/1"]) / length(a_rep)
		# compare expected and observed heterozygosities to calculate Fis
		fis <- (a_eh - a_het) / a_eh
		# add to output
    		output_rep <- c(popmap[a,1], "none", "fis", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), length(a_rep), length(a_rep[a_rep == "0/1"]), fis)
    		write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
	}
}


####################################
# helper function to determine frequency of 
# allele 1 (to be applied across rows)
# used in pi_tajima_theta function
####################################
p_function <- function(xxx) {
	xxxx <- c(substr(xxx, 1, 1), substr(xxx, 3, 3))
	xxxx <- length(xxxx[xxxx == 0])
	return(xxxx)
}


####################################
# helper function to determine frequency of 
# allele 2 (to be applied across rows)
# used in pi_tajima_theta function
####################################
q_function <- function(xxx) {
	xxxx <- c(substr(xxx, 1, 1), substr(xxx, 3, 3))
	xxxx <- length(xxxx[xxxx == 1])
	return(xxxx)
}


####################################
# calc watterson's theta, nucleotide 
# diversity, and tajima's d with the 
# formulae from Carlson et al. 2005: 
# 10.1101/gr.4326505 (equations 1, 
# 2, and 4) 
# calculated from VCF subsampled to
# a single population of interest
####################################
pi_tajima_theta <- function(xxx, popname, outname, chrom_info) {
	# remove sites with missing data
	for(a in 4:ncol(xxx)) {
		xxx <- xxx[xxx[,a] != "./.",]
	}
	# define number of chromosomes sampled
	n_chromosomes <- (ncol(xxx) - 3) * 2
	# define number of sites genotyped for this population without missing data
	n_sites <- nrow(xxx)
  
	# continue if there are sites genotyped
	if(n_sites > 0) {
		# subset to only polymorphic sites
		xxx <- xxx[apply(xxx[,4:ncol(xxx)], 1, polymorphic_function) > 1, ]
		
		# count polymorphic sites
		n_polymorphic <- nrow(xxx)
    
		# determine denominator of theta calculation
    		theta_denominator <- sum(1 / seq(from=1, to=(n_chromosomes - 1)))
    
		# continue if there are polymorphic sites
		if(n_polymorphic > 0) {
			# calculate per site theta
			theta <- (n_polymorphic / theta_denominator) / n_sites
      
			# calculate frequencies of alleles 1 and 2 for each polymorphic site
			p_freq <- as.vector(apply(xxx[,4:ncol(xxx)], 1, p_function)) / n_chromosomes
			q_freq <- as.vector(apply(xxx[,4:ncol(xxx)], 1, q_function)) / n_chromosomes
      
			# calculate nucleotide diversity (pi)
			pi <- sum(2 * p_freq * q_freq) * (n_chromosomes / (n_chromosomes - 1)) / n_sites
      
			# calculate D using the variance of d (and all necessary subcomponents)
			a1 <- sum(1 / seq(from=1, to=(n_chromosomes - 1)))
			a2 <- sum(1 / (seq(from=1, to=(n_chromosomes - 1))^2))
			b1 <- (n_chromosomes + 1) / (3 * (n_chromosomes - 1))
			b2 <- (2 * (n_chromosomes^2 + n_chromosomes + 3)) / (9 * n_chromosomes * (n_chromosomes - 1))
			c1 <- b1 - (1 / a1)
			c2 <- b2 - ((n_chromosomes + 2) / (a1 * n_chromosomes)) + (a2 / (a1^2))
			e1 <- c1 / a1
			e2 <- c2 / ((a1^2) + a2)
			Var_d <- (e1 * n_polymorphic) + (e2 * n_polymorphic * (n_polymorphic - 1))
			Tajima_D <- ((pi * n_sites - theta * n_sites)) / sqrt(Var_d)
      
			# write output
			output_rep1 <- c(popname, "none", "theta", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_polymorphic, theta)
    		write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
			output_rep2 <- c(popname, "none", "pi", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_polymorphic, pi)
    		write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
			output_rep3 <- c(popname, "none", "Tajima_D", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_polymorphic, Tajima_D)
    		write(output_rep3, file=outname, append=T, ncolumns=9, sep="\t")
   		} else { # if no polymorphic sites retained
   			output_rep1 <- c(popname, "none", "theta", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_polymorphic, "NA")
    		write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
			output_rep2 <- c(popname, "none", "pi", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_polymorphic, 0)
    		write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
			output_rep3 <- c(popname, "none", "Tajima_D", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_polymorphic, "NA")
    		write(output_rep3, file=outname, append=T, ncolumns=9, sep="\t")
   		} 
	}	else { # if no sites retained
		output_rep1 <- c(popname, "none", "theta", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), 0, 0, "NA")
    	write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
		output_rep2 <- c(popname, "none", "pi", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), 0, 0, 0)
    	write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
		output_rep3 <- c(popname, "none", "Tajima_D", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), 0, 0, "NA")
    	write(output_rep3, file=outname, append=T, ncolumns=9, sep="\t")
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
# helper function to identify 
# transition variants from biallelic 
# alleles in rows (used in titv function)
####################################
transitions <- function(x1) {
  a_rep <- sort(x1)
  if(a_rep[1] == "A" & a_rep[2] == "G") {
    return(1)
  } else if(a_rep[1] == "C" & a_rep[2] == "T") {
    return(1)
  } else {
    return(0)
  }
}


####################################
# calc transition / transversion 
# ratio from VCF
####################################
titv <- function(xxx, popmap, outname, chrom_info) {
	# filter to variable sites
	xxx <- xxx[xxx$ALT != ".", ]
	# sites used for titv 
	a_het <- nrow(xxx)
	# find number of transitions using helper function
	total_transitions <- sum(apply(xxx[,c(4,5)], 1, transitions))
	# calculate transition / transversion ratio
	titv_ratio <- total_transitions / (a_het - total_transitions)
	# output info
	output_rep1 <- c("all_inds", "none", "titv", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), a_het, a_het, titv_ratio)
    write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
}


####################################
# helper function to calc distance 
# matrix, from EEMS documentation
####################################
# Compute the diffs matrix using the "mean allele frequency"
# imputation method
bed2diffs_v2 <- function(genotypes) {
  
	nIndiv <- nrow(genotypes)
	nSites <- ncol(genotypes)
	missing <- is.na(genotypes)
  
	## Impute NAs with the column means (= twice the allele frequencies)
	geno_means <- colMeans(genotypes, na.rm = TRUE)
	# nIndiv rows of genotype means
	geno_means <- matrix(geno_means, nrow = nIndiv, ncol = nSites, byrow = TRUE) 
  
	## Set the means which correspond to observed genotypes to 0
	geno_means[missing == FALSE] <- 0
	## Set the missing genotypes to 0 (used to be NA) 
	genotypes[missing == TRUE] <- 0
	genotypes <- genotypes + geno_means
  
	similarities <- genotypes %*% t(genotypes) / nSites
	self_similarities <- diag(similarities)
	vector1s <- rep(1, nIndiv)
  
	diffs <- self_similarities %*% t(vector1s) + vector1s %*% t(self_similarities) - 2 * similarities
	diffs
}




####################################
# calc distance matrix
####################################
dist_mat <- function(xxx, outname) {
	# keep only genotypes
	xxx <- xxx[,10:ncol(xxx)]
	# convert genotypes to 0,1,2
	for(a in 1:ncol(xxx)) {
		xxx[xxx[,a] == "0/0",a] <- 0
		xxx[xxx[,a] == "0/1",a] <- 1
		xxx[xxx[,a] == "1/1",a] <- 2
		xxx[xxx[,a] == "./.",a] <- NA
		xxx[,a] <- as.numeric(xxx[,a])
	}
	# transpose for helper function
	xxx <- as.matrix(t(xxx))
	# calculate diff matrix
	diff_matrix <- bed2diffs_v2(xxx)
	n_snps <- ncol(xxx)
	outname <- paste0(outname, "__", n_snps, "__.diffs")
	# write output
	write.table(diff_matrix, file=outname, col.names=T, row.names=F, quote=F, sep="\t")
}


####################################
# calc private, shared, fixed 
# polymorphisms
####################################
poly_fps <- function(xxx1, xxx2, popname1, outname, chrom_info) {
	n_sites <- 0
	fixed_count <- 0
	shared_count <- 0
	private_count <- 0
	for(a in 1:nrow(xxx1)) {
		a_rep1 <- unique(as.character(unique(xxx1[a,])))
		a_rep2 <- unique(as.character(unique(xxx2[a,])))
		# remove missing
		a_rep1 <- a_rep1[a_rep1 != "./."]
		a_rep2 <- a_rep2[a_rep2 != "./."]
		a_rep <- unique(c(a_rep1, a_rep2))
		if(length(a_rep1) > 0 & length(a_rep) > 1) {
			if(length(a_rep1) == 1 & length(a_rep2) == 1) {
				fixed_count <- fixed_count + 1
				n_sites <- n_sites + 1
			} else if(length(a_rep1) == 1 & length(a_rep2) > 1) {
				# polymorphic outside focal population, don't count
			} else if(length(a_rep1) > 1 & length(a_rep2) == 1) {
				private_count <- private_count + 1
				n_sites <- n_sites + 1
			} else if(length(a_rep1) > 1 & length(a_rep2) > 1) {
				shared_count <- shared_count + 1
				n_sites <- n_sites + 1
			}
		}
	}
	output_rep <- c(popname1, "none", "private", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_sites, private_count)
    write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
	output_rep <- c(popname1, "none", "shared", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_sites, shared_count)
    write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
    output_rep <- c(popname1, "none", "fixed", chrom_info[1], as.numeric(chrom_info[2]),
    					as.numeric(chrom_info[3]), n_sites, n_sites, fixed_count)
    write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
}


####################################
# make FASTA from VCF
# heterozygous sites choose one allele
# at random
####################################
create_fasta_from_vcf <- function(xxx, outname, num_sites_needed) {
	# keep going only if enough sites
	if(num_sites_needed <= nrow(xxx)) {
		
		# subset the genotypes from the allele info
		allele_info <- xxx[,c(4,5)]
		genotypes <- xxx[,10:ncol(xxx)]
		# define names of individuals in output fasta
		output_names_fasta <- paste(">", colnames(genotypes), sep="")
		
		# define heterozygous ambiguities
    	het <- rep("N", nrow(allele_info))
    	het[allele_info[,1] == "A" & allele_info[,2] == "C"] <- "M"
     	het[allele_info[,1] == "A" & allele_info[,2] == "G"] <- "R"
      	het[allele_info[,1] == "A" & allele_info[,2] == "T"] <- "W"
       	het[allele_info[,1] == "C" & allele_info[,2] == "A"] <- "M"
       	het[allele_info[,1] == "C" & allele_info[,2] == "G"] <- "S"
       	het[allele_info[,1] == "C" & allele_info[,2] == "T"] <- "Y"
       	het[allele_info[,1] == "G" & allele_info[,2] == "A"] <- "R"
       	het[allele_info[,1] == "G" & allele_info[,2] == "C"] <- "S"
       	het[allele_info[,1] == "G" & allele_info[,2] == "T"] <- "K"
       	het[allele_info[,1] == "T" & allele_info[,2] == "A"] <- "W"
       	het[allele_info[,1] == "T" & allele_info[,2] == "C"] <- "Y"
       	het[allele_info[,1] == "T" & allele_info[,2] == "G"] <- "K"
       	
       	# convert all numbers in genotypes to actual bases and ambiguities		
		for(a in 1:ncol(genotypes)) {
			# keep only genotype
			genotypes[,a] <- substr(genotypes[,a], 1, 3)
			# convert to bases
			genotypes[,a][genotypes[,a] == "0/0"] <- allele_info[genotypes[,a] == "0/0",1]
			genotypes[,a][genotypes[,a] == "1/1"] <- allele_info[genotypes[,a] == "1/1",2]
			genotypes[,a][genotypes[,a] == "./."] <- "?"
			genotypes[,a][genotypes[,a] == "0/1"] <- het[genotypes[,a] == "0/1"]
			# choose random allele for heterozygotes
			genotypes[,a][genotypes[,a] == "R"] <- sample(c("A", "G"), length(genotypes[,a][genotypes[,a] == "R"]), replace=T)
			genotypes[,a][genotypes[,a] == "K"] <- sample(c("T", "G"), length(genotypes[,a][genotypes[,a] == "K"]), replace=T)
			genotypes[,a][genotypes[,a] == "S"] <- sample(c("C", "G"), length(genotypes[,a][genotypes[,a] == "S"]), replace=T)
			genotypes[,a][genotypes[,a] == "Y"] <- sample(c("C", "T"), length(genotypes[,a][genotypes[,a] == "Y"]), replace=T)
			genotypes[,a][genotypes[,a] == "M"] <- sample(c("A", "C"), length(genotypes[,a][genotypes[,a] == "M"]), replace=T)
			genotypes[,a][genotypes[,a] == "W"] <- sample(c("A", "T"), length(genotypes[,a][genotypes[,a] == "W"]), replace=T)

		}
    
		# write output
		for(a in 1:ncol(genotypes)) {
			if(a == 1) {
				write(output_names_fasta[a], file=outname, ncolumns=1)
				write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
			} else {
				write(output_names_fasta[a], file=outname, ncolumns=1, append=T)
				write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
			}
		}
	}
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
if(calc_heterozygosity == TRUE | calc_pi == TRUE | calc_diff == TRUE | calc_FIS == TRUE | calc_titv == TRUE | calc_polymorphisms == TRUE) {
	# define output name
	output_name <- paste(outname_base, "__stats.txt", sep="")
	# write output file
	write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=output_name, sep="\t")
}


####################################
# calc FIS
####################################
if(calc_FIS == TRUE) {
	# modify popmap to only have populations with N >= calc_FIS_minN
	pop_counts <- table(populations[,2])
	pop_counts <- names(pop_counts)[pop_counts >= calc_FIS_minN]
	populations2 <- populations[populations[,2] %in% pop_counts,]
	# calculate FIS per population
	fis(input_file, populations2, output_name, chromosome_info)
}


####################################
# calc heterozygosity
####################################
if(calc_heterozygosity == TRUE) {
	# calculate heterozygosity for each individual
	heterozygosity(input_file, output_name, chromosome_info)
}


####################################
# calc pi, theta, Tajima D
####################################
if(calc_pi == TRUE) {
	# modify popmap to only have populations with N >= calc_pi_minN
	pop_counts <- table(populations[,2])
	pop_counts <- names(pop_counts)[pop_counts >= calc_pi_minN]
	populations2 <- populations[populations[,2] %in% pop_counts,]
	
	# calculate theta, pi, and Tajima D for each population
	for(a in 1:length(unique(populations2[,2]))) {
		a_pop <- unique(populations2[,2])[a]
		# subset vcf input
		a_input <- cbind(input_file[,c(1,4,5)], 
			input_file[,colnames(input_file) %in% 
			populations2[populations2[,2] == a_pop, 1]])
		# calculate stats per pop
		pi_tajima_theta(a_input, a_pop, output_name, chromosome_info)
	}
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


####################################
# calc TI/TV
####################################
if(calc_titv == TRUE) {
	# calculate transition / transversion ratio
	titv(input_file, populations, output_name, chromosome_info)
}


####################################
# calc distance matrix
####################################
if(calc_genetic_distance == TRUE) {
	# calculate genetic distance matrix
	dist_mat(input_file, outname_base)
}


####################################
# calc fixed, private, shared
# polymorphisms
####################################
if(calc_polymorphisms == TRUE) {
	# modify popmap to only have populations with N >= calc_pi_minN
	pop_counts <- table(populations[,2])
	pop_counts <- names(pop_counts)[pop_counts >= calc_polymorphisms_minN]
	populations2 <- populations[populations[,2] %in% pop_counts,]
	
	# only polymorphic sites
	input_file2 <- input_file[input_file$ALT != ".", ]
	input_file2 <- input_file2[,10:ncol(input_file2)]
	
	# calculate polymorphism types for each population
	for(a in 1:length(unique(populations2[,2]))) {
		a_pop <- unique(populations2[,2])[a]
		# subset vcf input
		a_input1 <- input_file2[,colnames(input_file2) %in% 
			populations2[populations2[,2] == a_pop, 1]]
		a_input2 <- input_file2[,colnames(input_file2) %in% 
			populations2[populations2[,2] == a_pop, 1] == FALSE]
		# calculate per pop
		poly_fps(a_input1, a_input2, a_pop, output_name, chromosome_info)
	}
}


####################################
# make FASTA file
####################################
if(make_FASTA == TRUE) {
	# define output fasta name
	output_fasta <- paste(strsplit(args[1], ".recode.vcf")[[1]][1], ".fasta", sep="")
	# create fasta sequence alignments for all files with a minimum number of sites
	create_fasta_from_vcf(input_file, output_fasta, min_sites)
}







