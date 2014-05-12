# TODO: Add comment
# 
# Author: niko.popitsch@univie.ac.at
###############################################################################

library(ggplot2) 
library(reshape2)
library(gridExtra)
library(directlabels)
library(cwhmisc)


# input: 2 SNP datasets (calculated from VCF files)
data1=data.frame(read.csv("/project/oesi/genomicAmbiguity/snp-test/TEST-1T1N.union.vcf.csv", sep="\t", header=F,comment.char = "#"))
data2=data.frame(read.csv("/project/oesi/genomicAmbiguity/snp-test/TEST-2T2N.union.vcf.csv", sep="\t", header=F,comment.char = "#"))
colnames(data1)=c("qual", "score", "depth")
colnames(data2)=c("qual", "score", "depth")

# plot score vs. depth
g1=ggplot(data1, aes(score, depth)) + geom_point() + ggtitle("Exon SNP data set 1") + geom_smooth()
# plot score vs. quality
g2=ggplot(data1, aes(score, qual)) + geom_point() + geom_smooth()
# plost quality vs. depth
g3=ggplot(data1, aes(qual, depth)) + geom_point() + geom_smooth()
# plot score histogram
g4=ggplot(data1, aes(score)) + geom_histogram(binwidth=.1)

g5=ggplot(data2, aes(score, depth)) + geom_point() + ggtitle("Exon SNP data set 2")  + geom_smooth()
g6=ggplot(data2, aes(score, qual)) + geom_point() + geom_smooth()
g7=ggplot(data2, aes(qual, depth)) + geom_point() + geom_smooth()
g8=ggplot(data2, aes(score)) + geom_histogram(binwidth=.1)

pdf("/project/oesi/genomicAmbiguity/snp-test/TEST-correlations.pdf", width=10,height=16)
par(mfrow=c(1,1))
grid.arrange(g1, g5, g2, g6, g3, g7, g4, g8, nrow=4)
dev.off()