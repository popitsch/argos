
data1=read.table("c:/data/genomicAmbiguity/stackedHistPrec.csv", sep="\t", header=T, row.names=1, comment.char="#")

#rownames(data1)=c("<10" , "10-20" , " >20" )

#barplot(as.matrix(data1), main="ISS distribution", xlab="organism")

plot(data1$d.mel, log="y", type="l", col="blue", main="ISS distribution", xlab="ISS", ylab="log(%)" )
lines(data1$e.coli, col="red")
#lines(data1$h.sapiens, col="green")


