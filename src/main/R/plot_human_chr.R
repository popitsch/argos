#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomeGraphs")

# see https://www.biostars.org/p/378/

# prepareGenomePlot example
library(quantsmooth)
# construct genomic positions
CHR<-sample(22,40,replace=TRUE)  # Chromosomes
MapInfo<-lengthChromosome(CHR,"bases")*runif(length(CHR)) # position on chromosome
chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, organism="hsa")
# Chrompos returns a matrix with the positions of the elements on the plot
# You can use all kinds of base graphics functions to annotate the chromosomes
points(chrompos[,2],chrompos[,1]+0.1,pch="x",col="red")
# Show connection between 3rd and 4th element
segments(chrompos[3,2],chrompos[3,1],chrompos[4,2],chrompos[4,1],col="blue",lwd=2)