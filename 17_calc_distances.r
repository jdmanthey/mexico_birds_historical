
library(vcfR)
library(adegenet)
library(StAMPP)

# loop for each of the species to read in the vcf, create a distance matrix, and output in a
# format compatible with splitstree

x_files <- list.files(pattern="*structure.vcf")
for(a in 1:length(x_files)) {
	# read vcf
	a_rep <- read.vcfR(x_files[a])
	
	# convert to genlight
	a_rep <- vcfR2genlight(a_rep)
	
	# give fake population names (not used anyway)
	pop(a_rep) <- a_rep@ind.names
	
	# calculate Nei's distance 
	a_rep <- stamppNeisD(a_rep, pop = FALSE)
	
	new_names <- row.names(a_rep)
	new_names <- substr(new_names, nchar(new_names) - 7, nchar(new_names))
	new_names <- sapply(strsplit(new_names, "_"), "[[", 2)
	
	rownames(a_rep) <- new_names
	
	# output name
	outname <- paste0(strsplit(x_files[a], "\\.")[[1]][1], "_distmat.phy")
	
	# write output distance matrix in phylip format
	stamppPhylip(distance.mat=a_rep, file= outname)
	
}
