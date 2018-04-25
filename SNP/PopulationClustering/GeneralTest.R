#basic summary based on the genotyping of TCGA normal samples
#Date: April 16, 2018
setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/SNP/GeneralTest/")

################################Allele Frequency###################################
tcgafreq<-read.table("normal_sample_summary.frq",h=T)

hapmapfreq<-read.table("hapmap_summary.frq",h=T)

overlapsnp<-intersect(tcgafreq[,2],hapmapfreq[,2])

tcgamaf<-tcgafreq[pmatch(overlapsnp,tcgafreq[,2]),5]

png("MutationFrequencyDistribution.png",height=1024,width=1024,res=196)
plot(density(tcgamaf[-which(is.na(tcgamaf))]),xlab="Allele Frequency",main="TCGA-BRCA Allele Frequency")
dev.off()

hapmafmaf<-hapmapfreq[pmatch(overlapsnp,hapmapfreq[,2]),5]

png("tcgaVShapmap.png",height=1024,width=1024,res=196)
plot(tcgamaf[1:10000],hapmafmaf[1:10000],pch=19,cex=0.5,xlab="TCGA Maf (10000 SNP)",ylab="HapMap Maf (10000 SNP)")
lmt<- lm(hapmafmaf[1:10000]~tcgamaf[1:10000])#cor 0.84
abline(lmt,col="red")

################################Missing Value###################################
tcgaMissing<-read.table("normal_sample_summary.lmiss",h=T)

tcgamis<-tcgaMissing[,5]

freqC<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
misnumb<-sapply(freqC,function(x)length(which(tcgamis>x)))

MisMat<-as.matrix(rbind(freqC,misnumb,misnumb/905188))
write.csv(MisMat,"SNP_MissingRate.csv")

################################Missing Value###################################
hwefreq<-read.table("normal_sample_summary.hwe",h=T)
tcgahwef<-hwefreq[,9]

snphew<-hwefreq[which(tcgahwef<1e-06),2]
write.table(snphew,"HEW_failed_snps.txt",row.names=F,col.names=F,quote=F)

miscase<-read.table("missing_snps_across_cases.txt",h=F)
number<-sum(miscase[which(miscase[,2]/1150>0.1),1])#105 SNPs have a missing rate > 10% across 1150 samples

################################Sex Check###################################
sexpred<-read.table("normal_sample_summary.sexcheck.test",h=T)

tcgamate<-as.data.frame(read.csv("../../metadata/TCGA_BRCA_Clinical_Info.csv",h=T))#TCGA clinial info: RACE SAMPLE name SAMPLE ID
id<-read.table("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/TCGA_SNPSample_Annotation.txt",h=T)#sample name + TCGA ID

pcasplname<-gsub(".SNP6","",sexpred[,1])#sample name in PCA results

idsplname<-gsub(".birdseed.data.txt","",id[,2])#extract sample name in id

sampleID<-substr(id[pmatch(pcasplname,idsplname),7],1,12) #match sample name with TCGA Sample ID, etc."TCGA-E9-A24A"

sampleindex<-unlist(sapply(sampleID,function(x)if(length(which(tcgamate[,1]==x))!=0){return(which(tcgamate[,1]==x))}else{return(NA)}))

sex<-tcgamate$gender[sampleindex]

mismatchSample<-sexpred[-which(sexpred[,4]==c(2,1,0)[sex]),1]

sampleID<-substr(id[pmatch(gsub(".SNP6","",mismatchSample),idsplname),7],1,12) 




