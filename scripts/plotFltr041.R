rm(list=ls())

#setwd("C:/Users/cbird/Documents/GCL/scripts/fltrVCF/scripts")
library(tidyverse)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[3] = paste(args[1], ".pdf", sep="")
}

inputFile <- args[1]
lowPCT <- as.numeric(args[2])
highPCT <- as.numeric(args[3])
outputFile <- args[4]

#read in ldepth haps file
df <- read.table(inputFile, header=FALSE, sep='\t')
colnames(df) <- c("Contig", "BP", "Mean_Mean_Cvg")

cvgCutoffs <- quantile(df$Mean_Mean_Cvg, probs=c(lowPCT, highPCT))

maxCVG <- max(df$Mean_Mean_Cvg)

p1 <- ggplot(df, aes(x=Mean_Mean_Cvg)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("Contig Coverage") +
  scale_x_continuous(limits = c(0, maxCVG))

df1 <- df[df$Mean_Mean_Cvg >= cvgCutoffs[1] & df$Mean_Mean_Cvg <= cvgCutoffs[2],]

p2 <- ggplot(df1, aes(x=Mean_Mean_Cvg)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle(paste("Contig Coverage: ", round(cvgCutoffs[1],1), "< CVG <", round(cvgCutoffs[2],1), sep=" ")) +
  scale_x_continuous(limits = c(0, maxCVG))


# Output Plots
pdf(file=outputFile, width=8.5, height=11)
  grid.arrange(p1, p2, ncol=1, nrow=2)
dev.off()

