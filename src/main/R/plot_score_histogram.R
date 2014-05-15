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
#data1=data.frame(read.csv("c:/data/genomicAmbiguity/human/HIST-hg19-CCDS-ISS.csv", sep="\t", header=T,comment.char = "#"))

#par(mfrow=c(1,2))
#
#dataISSEx=data.frame(read.csv("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-EXONS-ISS.hist.csv", sep="\t", header=T,comment.char = "#"))
#colnames(dataISSEx)=c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100")
#barplot(as.matrix(dataISSEx), xlab="Avg. ISS score", ylab="frequency", main="Dros ISS exons")
#
#dataISSIn=data.frame(read.csv("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-INTRONS-ISS.hist.csv", sep="\t", header=T,comment.char = "#"))
#colnames(dataISSIn)=c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100")
#barplot(as.matrix(dataISSIn), xlab="Avg. ISS score", ylab="frequency", main="Dros ISS introns")
#
## NO need to plot AMB - the average is enough
#dataAMBEx=data.frame(read.csv("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-EXONS-AMB.csv", sep="\t", header=F,comment.char = "#"))
#print("average AMB score for all exons")
#print(ave(dataAMBEx[,6])[1])
#dataAMBIn=data.frame(read.csv("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-INTRONS-AMB.csv", sep="\t", header=F,comment.char = "#"))
#print("average AMB score for all introns")
#print(ave(dataAMBIn[,6])[1])

name="dros_partial/droso_r5_partial"
#name="hg19-NO_SCORE_LIMIT/hg19-NO_SCORE_LIMIT"

dataISSEx=data.frame(read.csv(paste("/project2/oesi/genAmb/output/",name,"-EXONS-ISS.csv",sep=""), sep="\t", header=F,comment.char = "#"))
dataISSIn=data.frame(read.csv(paste("/project2/oesi/genAmb/output/",name,"-INTRONS-ISS.csv",sep=""), sep="\t", header=F,comment.char = "#"))
dataISSMerged=data.frame(read.csv(paste("/project2/oesi/genAmb/output/",name,"-MERGED-ISS.csv",sep=""), sep="\t", header=F,comment.char = "#"))

#dataAMBEx=data.frame(read.csv("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-EXONS-AMB.csv", sep="\t", header=F,comment.char = "#"))
#dataAMBIn=data.frame(read.csv("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-INTRONS-AMB.csv", sep="\t", header=F,comment.char = "#"))


#pdf(paste("/project2/oesi/genAmb/output/",name,"-exonintron-density.pdf",sep=""), width=10,height=8)

par(mfrow=c(1,1))

#plot(density(dataISSEx$V6, bw=0.2),  col="blue")
#lines(density(dataISSIn$V6, bw=0.2),  col="red")
##abline(v=ave(dataISSEx$V6), col="blue")
##abline(v=ave(dataISSIn$V6), col="red")
#abline(v=median(dataISSEx$V6), col="blue", lty=2)
#abline(v=median(dataISSIn$V6), col="red", lty=2)
#
## show only the low scoring exons
#plot(density(dataISSEx[dataISSEx$V6<10,]$V6, bw=0.2),  col="blue")
#lines(density(dataISSIn[dataISSIn$V6<10,]$V6, bw=0.2),  col="red")
##abline(v=ave(dataISSEx[dataISSEx$V6<10,]$V6), col="blue")
##abline(v=ave(dataISSIn[dataISSIn$V6<10,]$V6), col="red")
#abline(v=median(dataISSEx[dataISSEx$V6<10,]$V6), col="blue", lty=2)
#abline(v=median(dataISSIn[dataISSIn$V6<10,]$V6), col="red", lty=2)

plot(dataISSMerged$V2, dataISSMerged$V3, xlim=c(0,100), ylim=c(0,100), col="blue", pch="o")
abline(0,1)
# list the genes with high exon ISS and low intron ISS
dataISSMerged[dataISSMerged$V2>80 & dataISSMerged$V3<20 & dataISSMerged$V3 >= 0,]

#dev.off()
