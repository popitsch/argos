#
# Plots a MAPQ analysis result
#
library(scales)
library(squash)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) < 1) {
	stop("usage: /software/R/R-current/bin/Rscript MapqAnalyzer.R [input-dat] [xlab] [ylab]")
}

fn=args[1]
xlab=args[2]
ylab=args[3]

#data=fread("/project2/oesi/genAmb-OLD/test_data/human_wgs/mapping",sep="\t", nrows=277076786, header=F)
raw=fread(fn, sep="\t", header=F)
dat = data.frame(raw)
output = paste( fn,"-",xlab,"-",ylab,".PLOT.pdf", sep="")

pdf(output, width=10,height=10)

#par(mfrow=c(1,2))
#plot(	dat$V2, 
#		dat$V3, 
#		pch=19, 
#		col = alpha("blue", 0.05),
#		xlab=xlab,
#		ylab=ylab
#		)

# draw headtmap		
hist2(dat$V2, 
		dat$V3, 
		xlab=xlab, 
		ylab=ylab,
		zsize=1,
		breaks=prettyLog,
		key.args = list(stretch = 3))		

#dat[dat[, "V2"]>60, ]
#dat[dat[, "V3"]>100, ]
#dat[dat[, "V2"]==60 & mapqvsiss[, "V3"]==0, ]

print('Finished.')
dev.off()