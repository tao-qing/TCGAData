#TCGA SNP Array Population Clustering Analysis based on Germline Mutation
#Date: April 13, 2018
setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/SNP/PopulationClustering/")

############################Population Clustering Analysis#################################
tcgamate<-as.data.frame(read.csv("../../metadata/TCGA_BRCA_Clinical_Info.csv",h=T))#TCGA clinial info: RACE SAMPLE name SAMPLE ID
id<-read.table("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/TCGA_SNPSample_Annotation.txt",h=T)#sample name + TCGA ID

pcamat<-read.table("merge_tcga_brca_2663_snp_normal_maf0.05_HEW_autosome_hapmapSNP.eigenvec")#PCA results
pcasplname<-gsub(".SNP6","",pcamat[,2])#sample name in PCA results

idsplname<-gsub(".birdseed.data.txt","",id[,2])#extract sample name in id

sampleID<-substr(id[pmatch(pcasplname,idsplname),7],1,12) #match sample name with TCGA Sample ID, etc."TCGA-E9-A24A"

sampleindex<-unlist(sapply(sampleID,function(x)which(tcgamate[,1]==x)))

race<-tcgamate$race_list[sampleindex]
col<-c("red","green","black","blue","purple")[as.factor(race)]
png("PCA_population_cluster_HapMapLocus.png",height = 1024,width = 1024,res=200)
plot(pcamat[,3],pcamat[,4],pch=19,col=col,cex=0.5,ylab="PC2",xlab="PC1")
legend("topright",c("White","Asian","Black","Indian","Unknown"),col=c("purple","black","blue","green","red"),pch=19)
dev.off()



############################HapMap only#################################
hapmappcamat<-read.table("hapmap_only.eigenvec")#PCA results
samplename<-hapmappcamat[,2]
hapmapinfo<-as.matrix(read.table("relationships_w_pops_051208.txt",sep="\t",h=T))
race<-hapmapinfo[pmatch(samplename,hapmapinfo[,2]),7]

col<-c("red","green","black","blue")[as.factor(race)]
png("Hapmap_PCA_population_cluster.png",height = 1024,width = 1024,res=200)
plot(hapmappcamat[,3],hapmappcamat[,4],pch=19,col=col,cex=0.5,ylab="PC2",xlab="PC1")
legend("topleft",c("ASW","CEU","MEX","YRI"),col=c("red","green","black","blue"),pch=19)
dev.off()

hapmapinfo<-read.table("hapmap.pop")

############################HCA Analysis#################################

hcamat<-read.table("merge_tcga_brca_2663_snp_cluster.cluster3",r=1)


############################Sample Origin Check##########################

forensic35Mat<-read.table("merge_tcga_brca_132forensicLocus.tped",h=T)
sampleName<-read.table("merge_tcga_brca_132forensicLocus.tfam")

rownames(forensic35Mat)<-forensic35Mat[,2]
forensic35Mat<-forensic35Mat[,-c(1:4)]
colnames(forensic35Mat)<-sampleName[,1]

SaplAnno<-read.csv("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/TCGA_Array6_Sample_information.csv")

samples<-unique(SaplAnno[,3])
frepMat<-sapply(samples,function(x){
        index<-grep(x,SaplAnno[,3])
        if(length(index)!=2){
          return(NA)
        }else{
          a<-forensic35Mat[,grep(SaplAnno[index[1],2],colnames(forensic35Mat))]
          b<-forensic35Mat[,grep(SaplAnno[index[2],2],colnames(forensic35Mat))]
          submat=cbind(a,b)
          freq<-length(which(apply(submat,1,function(x)x[1]==x[2])))/34
          return(freq)
        }
})


####
forensic35Mat<-read.table("merge_tcga_brca_132forensicLocus.ped",h=F,sep="\t")
sampleName<-read.table("merge_tcga_brca_132forensicLocus.tfam")

rownames(forensic35Mat)<-forensic35Mat[,1]
forensic35Mat<-forensic35Mat[,-c(1:6)]
#colnames(forensic35Mat)<-sampleName[,1]

SaplAnno<-read.csv("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/TCGA_Array6_Sample_information.csv")

samples<-unique(SaplAnno[,3])

#frepNumber变量为样本被检测的Array的数量。

frepNumber<-sapply(samples,function(x){
  index<-grep(x,SaplAnno[,3])
    return(length(index))#71个
})

frepMat<-sapply(samples,function(x){
  index<-grep(x,SaplAnno[,3])
  if(length(index)!=2){
    return(NA)#71个
  }else{
    a<-forensic35Mat[grep(SaplAnno[index[1],2],rownames(forensic35Mat)),]
    b<-forensic35Mat[grep(SaplAnno[index[2],2],rownames(forensic35Mat)),]
    submat=rbind(a,b)
    freq<-length(which(apply(submat,2,function(x)x[1]==x[2])))/34
    return(freq)
   }
})

samples[which(frepMat<0.6)]
frepMatrmna<-frepMat[-which(is.na(frepMat))]
percentage<-1-length(which(frepMatrmna>0.85))/1027
#about 5.8% of tumor and normal samples are not from the same individual 

library(ggplot2)

box_plot_mat<-as.data.frame(cbind(samples[-which(is.na(frepMat))],frepMatrmna))
colnames(box_plot_mat)<-c("Sample","Percentage")
ggplot(box_plot_mat, aes(x=Sample, y=Percentage, fill=Sample))+geom_boxplot()+guides(fill=FALSE)+ylab("Percentage")++theme(legend.position = "none",axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10)) 
ggsave("ForensicLocusOverlapPercentage.tiff",height=10,width=6,units="cm")

save.image("20180415.RData")

