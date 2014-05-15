# TODO: Add comment
# 
# Author: niko.popitsch@univie.ac.at
###############################################################################

# plot a diagram of the ReadScoresStatistics data

library(ggplot2) 
library(reshape2)
library(gridExtra)
library(directlabels)
library(cwhmisc)

data="/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19.scores.gz.hist.csv"
label="hg19, no score limit, local alignment"
		
# input: binned histogram
#data1=data.frame(read.csv("/project2/oesi/genAmb/output/ecK12/eck12_MG1655_ecoli-chr.scores.gz.hist.csv", sep="\t", header=T,comment.char = "#"))
data1=data.frame(read.csv(data, sep="\t", header=T,comment.char = "#"))
#colnames(data1)=c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100")

fileGraphOutput=paste(data,".pdf", sep="")

pdf(fileGraphOutput, width=10,height=8)

barplot(as.matrix(data1), xlab="# of scores", ylab="frequency", main=label)
maxval=max(data1)
maxbin=colnames(data1)[which.max(data1)]
text(1000,100000, paste("max:",maxval,"bin:",maxbin))

dev.off()

