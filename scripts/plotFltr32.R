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
hetCut <- as.numeric(args[2])
outputFile <- args[3]

#read in ldepth haps file
df <- read.table(inputFile, header=FALSE, sep='\t')
colnames(df) <- c("Contig", "SNPs", "Mean_Num_Het", "Mean_Num_Ind_Minus_Missing", "PropHet")

p1 <- ggplot(df, aes(x=PropHet)) +
  geom_histogram(color="black", fill="white", binwidth = 0.01) +
  theme_classic() +
  ggtitle("Mean Proportion Hetero Per Contig") +
  scale_x_continuous(limits = c(0, 1))

df1 <- df[df$PropHet <= hetCut,]

p2 <- ggplot(df1, aes(x=PropHet)) +
  geom_histogram(color="black", fill="white", binwidth = 0.01) +
  theme_classic() +
  ggtitle(paste("Mean Proportion Hetero Per Contig < ", hetCut, sep=" ")) +
  scale_x_continuous(limits = c(0, 1))

# Output Plots
pdf(file=outputFile, width=8.5, height=11)
  grid.arrange(p1, p2, ncol=1, nrow=2)
dev.off()

