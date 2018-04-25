#TCGA Sample Classification
#TCGA assign the tumor samples with 0-9, and normal samples with 10-19
setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/")

tcgamate<-as.data.frame(read.csv("./metadata/TCGA_BRCA_Clinical_Info.csv",h=T))#TCGA clinial info: RACE SAMPLE name SAMPLE ID
id<-read.table("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/TCGA_SNPSample_Annotation.txt",h=T)#sample name + TCGA ID
SaplAnno<-cbind(gsub(".birdseed.data.txt","",id[,2]),substr(id[,7],1,12),substr(id[,7],14,16),c("Tumor","Normal")[as.factor(as.numeric(substr(id[,7],14,15))>=10)])
write.csv(SaplAnno,"TCGA_Array6_Sample_information.csv")

tumorsample<-SaplAnno[which(SaplAnno[,4]=="Tumor"),]
write.table(paste0(tumorsample[,1],".SNP6"),"tumorsamplename.txt",quote=F,row.names=F,col.names=F)
normalsample<-SaplAnno[which(SaplAnno[,4]=="Normal"),]
write.table(paste0(normalsample[,1],".SNP6"),"normalsamplename.txt",quote=F,row.names=F,col.names=F)

#01A  01B  06A  10A  10B  10C  10D  11A  11B 