rm(list=ls())

#setwd("C:/Users/cbird/Documents/GCL/scripts/fltrVCF/scripts")
library(tidyverse)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[6] = paste(args[1], ".pdf", sep="")
}

inputFile <- args[1]
lowPCT_Mean_Mean_Cvg <- as.numeric(args[2])
highPCT_Mean_Mean_Cvg <- as.numeric(args[3])
lowPCT_Mean_CV_Cvg <- as.numeric(args[4])
highPCT_Mean_CV_Cvg <- as.numeric(args[5])
outputFile <- args[6]

#read in ldepth haps file
df <- read.table(inputFile, header=FALSE, sep='\t')
colnames(df) <- c("Contig", "BP", "Mean_Mean_Cvg", "CV_Mean_Cvg", "Mean_CV_Cvg", "CV_CV_Cvg")

Cutoffs_Mean_Mean_Cvg <- quantile(df$Mean_Mean_Cvg, probs=c(lowPCT_Mean_Mean_Cvg, highPCT_Mean_Mean_Cvg))
Cutoffs_Mean_CV_Cvg <- quantile(df$Mean_CV_Cvg, probs=c(lowPCT_Mean_CV_Cvg , highPCT_Mean_CV_Cvg ))

#note that I set thse to Mean_CV_Cvg because these do not exist in fltrVCF right
Cutoffs_CV_Mean_Cvg <- quantile(df$CV_Mean_Cvg, probs=c(lowPCT_Mean_CV_Cvg , highPCT_Mean_CV_Cvg ))
Cutoffs_CV_CV_Cvg <- quantile(df$CV_CV_Cvg, probs=c(lowPCT_Mean_CV_Cvg , highPCT_Mean_CV_Cvg ))

maxMean_Mean_Cvg <- max(df$Mean_Mean_Cvg)
maxMean_CV_Cvg <- max(df$Mean_CV_Cvg)
maxCV_Mean_Cvg <- max(df$CV_Mean_Cvg)
maxCV_CV_Cvg <- max(df$CV_CV_Cvg)

#make histograms for Mean_Mean_Cvg
p1_Mean_Mean_Cvg <- ggplot(df, aes(x=Mean_Mean_Cvg)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("Mean of Mean Cvg") +
  scale_x_continuous(limits = c(0, maxMean_Mean_Cvg))

df1 <- df[df$Mean_Mean_Cvg >= Cutoffs_Mean_Mean_Cvg[1] & df$Mean_Mean_Cvg <= Cutoffs_Mean_Mean_Cvg[2],]

p2_Mean_Mean_Cvg <- ggplot(df1, aes(x=Mean_Mean_Cvg)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle(paste("Mean of Mean Cvg: ", round(cvgCutoffs[1],1), "< CVG <", round(cvgCutoffs[2],1), sep=" ")) +
  scale_x_continuous(limits = c(0, maxMean_Mean_Cvg))

#make histograms for Mean_CV_Cvg
p1_Mean_CV_Cvg <- ggplot(df, aes(x=Mean_CV_Cvg)) +
  geom_histogram(color="black", fill="white", bins = 100) +
  theme_classic() +
  ggtitle("Mean of CV of Mean Cvg") +
  scale_x_continuous(limits = c(0, maxMean_CV_Cvg))

df2 <- df[df$Mean_CV_Cvg >= Cutoffs_Mean_CV_Cvg[1] & df$Mean_CV_Cvg <= Cutoffs_Mean_CV_Cvg[2],]

p2_Mean_CV_Cvg <- ggplot(df2, aes(x=Mean_CV_Cvg)) +
  geom_histogram(color="black", fill="white", bins = 100) +
  theme_classic() +
  ggtitle(paste("Mean of CV of Mean Cvg: 0 < CVG <", round(Cutoffs_Mean_CV_Cvg[2],1), sep=" ")) +
  scale_x_continuous(limits = c(0, maxMean_CV_Cvg))

#make histograms for CV_Mean_Cvg
p1_CV_Mean_Cvg <- ggplot(df, aes(x=CV_Mean_Cvg)) +
  geom_histogram(color="black", fill="white", bins = 100) +
  theme_classic() +
  ggtitle("CV of Mean of Mean Cvg, Currently Not Filterable") +
  scale_x_continuous(limits = c(0, maxCV_Mean_Cvg))

df3 <- df[df$CV_Mean_Cvg >= Cutoffs_CV_Mean_Cvg[1] & df$CV_Mean_Cvg <= Cutoffs_CV_Mean_Cvg[2],]

p2_CV_Mean_Cvg <- ggplot(df3, aes(x=CV_Mean_Cvg)) +
  geom_histogram(color="black", fill="white", bins = 100) +
  theme_classic() +
  ggtitle(paste("Currently Not Filterable: CV of CV of Mean Cvg: 0 < CVG <", round(Cutoffs_CV_Mean_Cvg[2],1), sep=" ")) +
  scale_x_continuous(limits = c(0, maxCV_Mean_Cvg))

#make histograms for Mean_CV_Cvg
p1_CV_CV_Cvg <- ggplot(df, aes(x=CV_CV_Cvg)) +
  geom_histogram(color="black", fill="white", bins = 100) +
  theme_classic() +
  ggtitle("CV of CV of Mean Cvg, Currently Not Filterable") +
  scale_x_continuous(limits = c(0, maxCV_CV_Cvg))

df4 <- df[df$CV_CV_Cvg >= Cutoffs_CV_CV_Cvg[1] & df$CV_CV_Cvg <= Cutoffs_CV_CV_Cvg[2],]

p2_CV_CV_Cvg <- ggplot(df4, aes(x=CV_CV_Cvg)) +
  geom_histogram(color="black", fill="white", bins = 100) +
  theme_classic() +
  ggtitle(paste("Currently Not Filterable: CV of CV of Mean Cvg: 0 < CVG <", round(Cutoffs_CV_CV_Cvg[2],1), sep=" ")) +
  scale_x_continuous(limits = c(0, maxCV_CV_Cvg))


# Output Plots
pdf(file=outputFile, width=8.5, height=11)
  grid.arrange(p1_Mean_Mean_Cvg, p2_Mean_Mean_Cvg, ncol=1, nrow=2)
  grid.arrange(p1_Mean_CV_Cvg, p2_Mean_CV_Cvg, ncol=1, nrow=2)
  grid.arrange(p1_CV_Mean_Cvg, p2_CV_Mean_Cvg, ncol=1, nrow=2)
  grid.arrange(p1_CV_CV_Cvg, p2_CV_CV_Cvg, ncol=1, nrow=2)
dev.off()

