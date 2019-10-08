#!/usr/bin/env Rscript
#script to create graphs from fltrVCFstats.sbatch

#to run:

#Rscript plot_fltrVCFstats.R inputfile outputfile

rm(list=ls())

#setwd("C:/Users/cbird/Documents/GCL/scripts/fltrVCF/scripts")
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste(args[1], ".pdf", sep="")
}

inputFile <- args[1]
outputFile <- args[2]

#read in data
df <- read.table(inputFile, header=TRUE, sep='\t')
df$File <- factor(df$File, levels = as.character(df$File))
df$FilterID <- factor(df$FilterID, levels = as.character(df$FilterID))
df$totalGenotypes <- rowSums(df[7:12])
df$PropMissingGeno <- df$NumMissingGeno / df$totalGenotypes
df$PropGenoLess10X <- df$NumGenoLess10X / df$totalGenotypes
df$PropGeno10.19X <- df$NumGeno10.19X  / df$totalGenotypes
df$PropGeno20.49X <- df$NumGeno20.49X  / df$totalGenotypes
df$PropGeno50.99X <- df$NumGeno50.99X  / df$totalGenotypes
df$PropGeno100.999X <- df$NumGeno100.999X  / df$totalGenotypes

#stack data
df_cvgNum <- data.frame(df[1:3], stack(df[7:12]))
colnames(df_cvgNum) <- c("File", "FilterID", "FilterOrder", "NumGenotypes", "Depth")

df_cvgProp <- data.frame(df[1:3], stack(df[14:19]))
colnames(df_cvgProp) <- c("File", "FilterID", "FilterOrder", "ProportionGenotypes", "Depth")


#plot function
lineplot <- function(DATA=df, X='FilterID', Y='NumInd', G=1){
  ggplot(data=df, aes_string(x=X, y=Y, group=G)) +
  geom_line(color="black") +
  geom_point() +
  theme_classic()
}

#Create Plots
plotNumInd <- lineplot()
plotNumContigs <- lineplot(Y='NumContigs')
plotNumSNPs <- lineplot(Y='NumSNPs')

plotNumGenoDepth <- ggplot(data=df_cvgNum, aes(x=FilterID, y=NumGenotypes, color=Depth, group=Depth)) +
  geom_line() +
  geom_point() +
  theme_classic()

plotPropGenoDepth <- ggplot(data=df_cvgProp, aes(x=FilterID, y=ProportionGenotypes, color=Depth, group=Depth)) +
  geom_line() +
  geom_point() +
  theme_classic()

pdf(outputFile)
plotNumInd
plotNumContigs
plotNumSNPs
plotNumGenoDepth
plotPropGenoDepth
dev.off()
