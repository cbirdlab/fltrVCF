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
BPcutoff <- as.numeric(args[2])
outputFile <- args[3]

#read in ldepth file
#df <- read.table(inputFile, header=FALSE, sep='\t')
#colnames(df) <- c("Contig", "Position", "Mean_Cvg", "Var_Cvg")

#read in ldepth haps file
df <- read.table(inputFile, header=FALSE, sep='\t')
colnames(df) <- c("Contig", "BP", "Mean_Mean_Cvg")

maxBP <- max(df$BP)
minBP <- min(df$BP)

#plot function
scatplot <- function(DATA=df, X='BP', Y='Mean_Mean_Cvg', G=1, T=paste("Before Filter 31, Remove Contigs Shorter than", BPcutoff, "BP", sep=" "), YL="Contig Means of Mean Cvg Per Site"){
  ggplot(data=DATA, aes_string(x=X, y=Y, group=G)) +
    geom_point() +
    geom_smooth() +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    ggtitle(T) +
    ylab(YL) +
    xlab("Contig Length (BP)")
}

boxplots <- function(DATA=df, X='BP', Y='Mean_Mean_Cvg', T=paste("Before Filter 31, Remove Contigs Shorter than", BPcutoff, "BP", sep=" "), YL="Contig Means of Mean Cvg Per Site"){
  ggplot(data=DATA, aes_string(x=X, y=Y, group=X)) +
    geom_boxplot() +
    #geom_smooth() +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(T) +
    ylab(YL) +
    xlab("Contig Length (BP)")
}

#Create Hypothetical Filter Plots
scatList <- list()
boxList <- list()
hypCutoffs <- seq(minBP, maxBP, round((maxBP-minBP)/5, 0))

for(i in hypCutoffs){
  df1 <- df[df$BP > i,]
  #print(scatplot(DATA=df1, T=paste("Hypothetical Filter: Contigs Shorter Than", i, "BP are Removed", sep=" ")))
  #print(boxplots(DATA=df1, T=paste("Hypothetical Filter: Contigs Shorter Than", i, "BP are Removed", sep=" ")))
  p1 <- scatplot(DATA=df1, T=paste("Hypothetical Filter:", i, "BP", sep=" "), YL="Mean Cvg")
  p2 <- boxplots(DATA=df1, T=paste("Hypothetical Filter:", i, "BP", sep=" "), YL="Mean Cvg")
  scatList[[which(hypCutoffs == i)]] <- p1
  boxList[[which(hypCutoffs == i)]] <- p2
}


#Create Filter Plots
df1 <- df[df$BP > BPcutoff,]
scat1 <- scatplot(DATA=df1, T=paste("After Filter 31. Removed Contigs Shorter Than", BPcutoff, "BP", sep=" "))
box1 <- boxplots(DATA=df1, T=paste("After Filter 31. Removed Contigs Shorter Than", BPcutoff, "BP", sep=" "))

# Output Plots
pdf(file=outputFile, width=8.5, height=11)
  #actual filter plot
  grid.arrange(scat1, box1, ncol=1, nrow=2)
  #hypothetical filter plots
  grid.arrange(scatList[[1]], scatList[[2]],scatList[[3]], scatList[[4]], scatList[[5]], ncol=1, nrow=5)
  grid.arrange(boxList[[1]], boxList[[2]], boxList[[3]], boxList[[4]], boxList[[5]], ncol=1, nrow=5)
  
dev.off()

