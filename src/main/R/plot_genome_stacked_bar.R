
data1=read.table("/project2/oesi/genAmb/output/POSTER_DATA/stackedHistPrec.csv", sep="\t", header=T, row.names=1, comment.char="#")

#rownames(data1)=c("<10" , "10-20" , " >20" )

#barplot(as.matrix(data1), main="ISS distribution", xlab="organism")

plot(data1$e.coli, type="l", col="blue", main="ISS distribution", xlab="ISS", ylab="%" )
lines(data1$d.mel, col="red")
#lines(data1$h.sapiens, col="green")


