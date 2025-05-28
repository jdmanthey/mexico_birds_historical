options(scipen=999)

x <- scan("genotype_stats.txt", what="character", sep="\t")

n_inds <- 208

id <- x[1:n_inds * 7 - 6]
sites6 <- x[1:n_inds * 7 - 4]
sites8 <- x[1:n_inds * 7 - 2]
sites10 <- x[1:n_inds * 7 ]

output <- data.frame(id=as.character(id), sites6=as.numeric(sites6), sites8=as.numeric(sites8), sites10=as.numeric(sites10))

atl <- output[1:72,]
car <- output[73:106,]
cat <- output[107:162,]
myi <- output[163:208,]
atl <- atl[order(atl$sites6, decreasing=T),]
car <- car[order(car$sites6, decreasing=T),]
cat <- cat[order(cat$sites6, decreasing=T),]
myi <- myi[order(myi$sites6, decreasing=T),]

par(mfrow=c(1,2))
plot(output[,2:3], xlab="Sites Genotyped (Min. Cov. = 6)", ylab="Sites Genotyped (Min. Cov. = 8)")
plot(output[,c(2,4)], xlab="Sites Genotyped (Min. Cov. = 6)",  ylab="Sites Genotyped (Min. Cov. = 10)")

par(mfrow=c(1,4))
plot(atl[,2], ylim=c(0,1100000000), ylab="Atlapetes Sites Genotyped", xlab="Individuals Sorted by Sites Genotyped")
plot(car[,2], ylim=c(0,1100000000), ylab="Cardellina Sites Genotyped", xlab="Individuals Sorted by Sites Genotyped")
plot(cat[,2], ylim=c(0,1100000000), ylab="Catharus Sites Genotyped", xlab="Individuals Sorted by Sites Genotyped")
plot(myi[,2], ylim=c(0,1100000000), ylab="Myioborus Sites Genotyped", xlab="Individuals Sorted by Sites Genotyped")

