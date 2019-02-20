#!/usr/bin/Rscript
#this will make histograms of vcftools -geno-r2 output, which is for linkage disequilibrium
ld <- read.table("out.geno.ld", header=TRUE, sep="\t")
NAs <- ld == "-nan"
ld[NAs] <- NA
linked <- ld[which(ld[,5] > 0.1),]

pdf(file="ldhist.pdf")
hist(ld[,5])
hist(ld[which(ld[,5] > 0.1),5])
hist(ld[which(ld[,5] > 0.2),5])
hist(ld[which(ld[,5] > 0.3),5])
hist(ld[which(ld[,5] > 0.4),5])
hist(ld[which(ld[,5] > 0.5),5])
heatmap()
dev.off()
ld[which(ld[,5]>0.9),]

