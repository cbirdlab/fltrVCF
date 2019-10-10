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
snpFile <- args[2]  #"PIRE_SiganusSpinus.Z.15.15.Fltr31.1.snp" 
insFile <- args[3]  #"PIRE_SiganusSpinus.Z.15.15.Fltr31.1.ins" 
delFile <- args[4]  #"PIRE_SiganusSpinus.Z.15.15.Fltr31.1.del" 
lowCutoff <- as.numeric(args[5])
highCutoff <- as.numeric(args[6])
outputFile <- args[7]

#read in ldepth file
df <- read.table(inputFile, header=FALSE, sep='\t')
colnames(df) <- c("Contig", "Pos", "Mean_Cvg")
df$Pos <- as.integer(df$Pos)

dfsnp <- read.table(snpFile, header=FALSE, sep='\t')
colnames(dfsnp) <- c("Contig", "Pos", "ID", "Ref")
dfsnp$Pos <- as.integer(dfsnp$Pos)

dfins <- read.table(insFile, header=FALSE, sep='\t')
colnames(dfins) <- c("Contig", "Pos", "ID", "Ref")
dfins$Pos <- as.integer(dfins$Pos)

dfdel <- read.table(delFile, header=FALSE, sep='\t')
colnames(dfdel) <- c("Contig", "Pos", "ID", "Ref")
dfdel$Pos <- as.integer(dfdel$Pos)

maxPos <- max(df$Pos)
minPos <- min(df$Pos)

snp1 <- ggplot(dfsnp, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("SNPs")
snp2 <- ggplot(dfsnp, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("SNPs") +
  scale_x_continuous(limits = c(0, 50))
snp3 <- ggplot(dfsnp, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("SNPs") +
  facet_grid(vars(Ref))
ins1 <- ggplot(dfins, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("Insertions")
ins2 <- ggplot(dfins, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("Insertions") +
  scale_x_continuous(limits = c(0, 50))
del1 <- ggplot(dfdel, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("Deletions")
del2 <- ggplot(dfdel, aes(x=Pos)) +
  geom_histogram(color="black", fill="white", binwidth = 1) +
  theme_classic() +
  ggtitle("Deletions") +
  scale_x_continuous(limits = c(0, 50))


#plot function
boxplots <- function(DATA=df, X='Pos', Y='Mean_Cvg', T="Coverage by Position", YL="Mean Cvg Per Site", X1=minPos, X2=maxPos){
  ggplot(data=DATA, aes_string(x=X, y=Y, group=X)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(T) +
    ylab(YL) +
    xlab("Position") +
    scale_x_continuous(limits = c(X1, X2))
}

panels <- seq(1,ceiling(maxPos/150))
cvgList <- list()
for(i in panels){
  cvg1 <- boxplots(X1=0+(150*(i-1)), X2=150*i)
  cvgList[[i]] <- cvg1
}
cvg2 <- boxplots(X1=minPos, X2=minPos+50)
cvg3 <- boxplots(X1=maxPos-50, X2=maxPos)


# Output Plots

pdf(file=outputFile, width=8.5, height=11)
  if (length(panels) == 1) {
    print(cvgList[[1]])
  } else {
    grid.arrange(cvgList[[1]], cvgList[[2]] , ncol=1, nrow=2)
    if (length(panels) == 3) {
      print(cvgList[[3]])
    } else if (length(panels) > 3) {
      grid.arrange(cvgList[[3]], cvgList[[4]] , ncol=1, nrow=2)
      if (length(panels) == 5) {
        print(cvgList[[5]])
      } else if (length(panels) > 5) {
        grid.arrange(cvgList[[5]], cvgList[[6]] , ncol=1, nrow=2)
        if (length(panels) == 7) {
          print(cvgList[[7]])
        } else if (length(panels) > 7) {
          grid.arrange(cvgList[[7]], cvgList[[8]] , ncol=1, nrow=2)
          if (length(panels) == 9) {
            print(cvgList[[9]])
          } else if (length(panels) > 9) {
            grid.arrange(cvgList[[9]], cvgList[[10]] , ncol=1, nrow=2)
            if (length(panels) == 11) {
              print(cvgList[[11]])
            } else if (length(panels) > 11) {
              grid.arrange(cvgList[[11]], cvgList[[12]] , ncol=1, nrow=2)
              if (length(panels) == 13) {
                print(cvgList[[13]])
              } else if (length(panels) > 11) {
                print("Modify script, sequences too long")
              }
            }
          }
        }
      }
    }
    
  }
  grid.arrange(cvg2, cvg3, ncol=1, nrow=2)
  grid.arrange(snp1, snp2 , ncol=1, nrow=2)
  print(snp3)
  grid.arrange(ins1, ins2, ncol=1, nrow=2)
  grid.arrange(del1, del2, ncol=1, nrow=2)
dev.off()

