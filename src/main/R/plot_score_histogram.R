# TODO: Add comment
# 
# Author: niko.popitsch@univie.ac.at
###############################################################################

library(ggplot2) 
library(reshape2)
library(gridExtra)
library(directlabels)
library(cwhmisc)


# input: binned histogram
data1=data.frame(read.csv("c:/data/genomicAmbiguity/human/HIST-hg19-CCDS-ISS.csv", sep="\t", header=T,comment.char = "#"))
colnames(data1)=c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100")

barplot(as.matrix(data1), xlab="avg. score", ylab="frequency")



