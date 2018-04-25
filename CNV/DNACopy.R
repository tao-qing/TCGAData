
setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/CNV")
library(data.table)

mat<-as.data.frame(fread("test_tangent"))
matrmna<-mat[-which(is.na(mat$test_tangent_1)),]
matrmna$test_tangent_1<-round(log2((matrmna$test_tangent_1)/2),digits=4)

library(DNAcopy)


CNA.object <- CNA(cbind(matrmna$test_tangent_1),matrmna$CHROM,matrmna$POS,data.type="logratio",sampleid="test_tangent_1")

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)


png("cnv_test.png",width=2000,height=500)
plot(segment.smoothed.CNA.object, plot.type="w")
dev.off()



library(DNAcopy)
data(coriell)

