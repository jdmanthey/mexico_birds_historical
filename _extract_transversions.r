options(scipen=999)

# read in a small vcf (don't use for large vcf files)
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}

args <- commandArgs(trailingOnly = TRUE)

# read in vcf  
x <- read_vcf(args[1])

# output name
output_name <- args[2]

# remove transitions
x <- x[(x$REF == "C" & x$ALT == "T") == FALSE,]
x <- x[(x$REF == "T" & x$ALT == "C") == FALSE,]
x <- x[(x$REF == "G" & x$ALT == "A") == FALSE,]
x <- x[(x$REF == "A" & x$ALT == "G") == FALSE,]

# write to output that already has header
write.table(x, file=output_name, quote=F, col.names=F, row.names=F, sep="\t", append=T)



