rm(list=ls())

library(ggplot2)
library(gridExtra)

#setwd("C:/Users/cbird/Downloads")


args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste(args[1], ".pdf", sep="")
}

inputFile <- args[1] # "OpihiSK2017_EastMaui.C.Plate2.5.5.Fltr19.ind_stats.out"
outputFile <- args[2] 

df <- read.table(inputFile, header=TRUE)

#text to columns by delimiter
df2 <- data.frame(do.call('rbind',strsplit(as.character(df$Ind), '_', fixed=TRUE)))
colnames(df2) <- c("SampleGroup", "IndividualID")
# df3 <- data.frame(do.call('rbind',strsplit(as.character(df2$SampleGroup), '-', fixed=TRUE)))
# colnames(df3) <- c("Proj", "Region", "Site", "LifeStage")
df4 <- data.frame(do.call('rbind',strsplit(as.character(df2$IndividualID), '-', fixed=TRUE)))
#colnames(df4) <- c("IndID", "NGS_Library", "Lanes")
colnames(df4) <- c("IndID", "NGS_Library")

df5 <- cbind(df, df2, df4)


pdf(outputFile)

# histogram of possible paralogs
hp1 <- ggplot(data=df5, aes(x=Poss_Paralogs, fill=SampleGroup)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Histogram of Paralogous Contigs Per Individual") +
  xlab("Number of Contigs Identified as Poss_Paralogs Per Individual") +
  ylab("Number of Individuals")
hp2 <- ggplot(data=df5, aes(x=Poss_Paralogs, fill=NGS_Library)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Histogram of Paralogous Contigs Per Individual") +
  xlab("Number of Contigs Identified as Poss_Paralogs Per Individual") +
  ylab("Number of Individuals")

hl1 <- ggplot(data=df5, aes(x=Low_Coverage.Errors, fill=SampleGroup)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Histogram of LowCvg or Erroneous Contigs Per Individual") +
  xlab("Number of Contigs Identified as LowCvg or Erroneous Per Individual") +
  ylab("Number of Individuals")
hl2 <- ggplot(data=df5, aes(x=Low_Coverage.Errors, fill=NGS_Library)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Histogram of Paralogous Contigs Per Individual") +
  xlab("Number of Contigs Identified as LowCvg or Erroneous Per Individual") +
  ylab("Number of Individuals")

hm1 <- ggplot(data=df5, aes(x=Miss_Genotype, fill=SampleGroup)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Histogram of Contigs With Missing Genotypes Per Individual") +
  xlab("Number of Contigs With Missing Genotypes Per Individual") +
  ylab("Number of Individuals")
hm2 <- ggplot(data=df5, aes(x=Miss_Genotype, fill=NGS_Library)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Histogram of Contigs With Missing Genotypes Per Individual") +
  xlab("Number of Contigs With Missing Genotypes Per Individual") +
  ylab("Number of Individuals")

grid.arrange(hp1, hp2, ncol=1, nrow=2)
grid.arrange(hl1, hl2, ncol=1, nrow=2)
grid.arrange(hm1, hm2, ncol=1, nrow=2)
#grid.arrange(hp1, hp2, hl1, hl2, hm1, hm2, ncol=2, nrow=3)

# plot Paralogs by SampleGroup
bp1 <- ggplot(data=df5, aes(x=SampleGroup, y=Poss_Paralogs)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplots of Paralogous Contigs Per Individual by Sample Group") +
  xlab("Sample Group") +
  ylab("Number of Poss_Paralogs per Individual")
bp2 <- ggplot(data=df5, aes(x=SampleGroup, y=Poss_Paralogs, color=NGS_Library)) +
  geom_point(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Plot of Paralogous Contigs Per Individual by Sample Group & NGS Library") +
  xlab("Sample Group") +
  ylab("Number of Poss_Paralogs per Individual")

bl1 <- ggplot(data=df5, aes(x=SampleGroup, y=Low_Coverage.Errors)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplots of LowCvg or Erroneous Contigs Per Individual by Sample Group") +
  xlab("Sample Group") +
  ylab("Number of LowCvg or Erroneous Contigs per Individual")
bl2 <- ggplot(data=df5, aes(x=SampleGroup, y=Low_Coverage.Errors, color=NGS_Library)) +
  geom_point(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Plot of LowCvg or Erroneous Contigs Per Individual by Sample Group & NGS Library") +
  xlab("Sample Group") +
  ylab("Number of LowCvg or Erroneous Contigs per Individual")

bm1 <- ggplot(data=df5, aes(x=SampleGroup, y=Miss_Genotype)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplots of Contigs With Missing Genotypes Per Individual by Sample Group") +
  xlab("Sample Group") +
  ylab("Number of Contigs With Missing Genotypes per Individual")
bm2 <- ggplot(data=df5, aes(x=SampleGroup, y=Miss_Genotype, color=NGS_Library)) +
  geom_point(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Plot of Contigs With Missing Genotypes Per Individual by Sample Group & NGS Library") +
  xlab("Sample Group") +
  ylab("Number of Contigs With Missing Genotypes per Individual")

#grid.arrange(bp1, bp2, ncol=1, nrow=2)
#grid.arrange(bl1, bl2, ncol=1, nrow=2)
#grid.arrange(bm1, bm2, ncol=1, nrow=2)
grid.arrange(bp1,bl1,bm1, ncol=1, nrow=3)
#grid.arrange(bp2,bl2,bm2, ncol=1, nrow=3)

# plot Paralogs by PlatePool identity (lab effect)
bp3 <- ggplot(data=df5, aes(x=NGS_Library, y=Poss_Paralogs)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplots of Paralogous Contigs Per Individual by NGS Library") +
  xlab("NGS Library") +
  ylab("Number of Poss_Paralogs per Individual")
bp4 <- ggplot(data=df5, aes(x=NGS_Library, y=Poss_Paralogs, color=SampleGroup)) +
  geom_point(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Plot of Paralogous Contigs Per Individual by NGS Library & Sample Group") +
  xlab("NGS Library") +
  ylab("Number of Poss_Paralogs per Individual")

bl3 <- ggplot(data=df5, aes(x=NGS_Library, y=Low_Coverage.Errors)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplots of LowCvg or Erroneous Contigs Per Individual by NGS Library") +
  xlab("NGS Library") +
  ylab("Number of LowCvg or Erroneous Contigs per Individual")
bl4 <- ggplot(data=df5, aes(x=NGS_Library, y=Low_Coverage.Errors, color=SampleGroup)) +
  geom_point(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Plot of LowCvg or Erroneous Contigs Per Individual by NGS Library & Sample Group") +
  xlab("NGS Library") +
  ylab("Number of LowCvg or Erroneous Contigs per Individual")

bm3 <- ggplot(data=df5, aes(x=NGS_Library, y=Miss_Genotype)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplots of Contigs With Missing Genotypes Per Individual by NGS Library") +
  xlab("NGS Library") +
  ylab("Number of  Contigs With Missing Genotypes per Individual")
bm4 <- ggplot(data=df5, aes(x=NGS_Library, y=Miss_Genotype, color=SampleGroup)) +
  geom_point(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Plot of  Contigs With Missing Genotypes Per Individual by NGS Library & Sample Group") +
  xlab("NGS Library") +
  ylab("Number of  Contigs With Missing Genotypes per Individual")

grid.arrange(bp3,bl3,bm3, ncol=1, nrow=3)
#grid.arrange(bp4,bl4,bm4, ncol=1, nrow=3)

ggplot(data=df5, aes(x=Low_Coverage.Errors, y=Poss_Paralogs, color=SampleGroup)) +
  geom_point(size=5) +
  geom_smooth(aes(group=1), color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Scatterplot of (Paralogous Contigs Per Individual) vs (Contigs with LowCvg or Genotyping Errors Per Individual)") +
  xlab("Number of Low Cvg or Erroneous Contigs Per Individual") +
  ylab("Number of Poss_Paralogs per Individual") 


ggplot(data=df5, aes(x=Miss_Genotype, y=Poss_Paralogs, color=SampleGroup)) +
  geom_point(size=5) +
  geom_smooth(aes(group=1), color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Scatterplot of (Paralogous Contigs Per Individual) vs (Contigs with LowCvg or Genotyping Errors Per Individual)") +
  xlab("Number of Contigs w/ Missing Genotypes Per Individual") +
  ylab("Number of Poss_Paralogs per Individual")

dev.off()

