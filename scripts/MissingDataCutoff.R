#Before you can import your out.miss file to R, you must change it to a .csv file.
#This can be done by opening your out.imiss file in Excel and saving as .csv. 

#Load the ggplot2 package, and set your working directory to the location of your
#new .csv file. 
setwd(dir = "Downloads")
#Read in your .csv file as a data frame, substituing "Your_out.miss.csv" with the 
#your file name. 
OUT.MISS <- read.csv("Your_out.miss.csv", header=TRUE)
#Sort individuals by missing data. 
OUT.MISS$INDV <- factor(OUT.MISS$INDV, levels = OUT.MISS$INDV[order(OUT.MISS$F_MISS)])
#Make a generic plot of your missing data and select your cutoff. 
ggplot(OUT.MISS, aes(x=INDV, y=F_MISS)) + geom_point()
#Generate horizontal line data with your final cutoff. Change the yint and lt values
#to your selected cutoff for missing data. 
hline <- data.frame(yint=0.15, lt="0.15")
#Plot your missing data with your cutoff. 
ggplot(OUT.MISS, aes(x=INDV, y=F_MISS)) + geom_point() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, size=0.5, colour = "black")) + ggtitle("Individuals Sorted By Missing Data") + xlab("Individuals") + ylab("Frequency of Missing Data") +geom_hline(data =hline, aes(yintercept =yint, linetype=lt), color="red", show.legend=TRUE) + scale_colour_discrete(guide="none") + scale_linetype_manual(name = 'Cutoff', values=1, guide="legend")
