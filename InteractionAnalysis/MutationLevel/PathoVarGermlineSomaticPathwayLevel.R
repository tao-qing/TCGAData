setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/InteractionAnalysis/MutationLevel")

library(data.table)
############################################Breast Cancer######################################################
germSNV<-read.csv("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/Data/PCAGermlineMutation/PCA_pathVar_mutations_BRCA.maf",sep="\t")
somaticSNV<-read.csv("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/Data/PCATCGAsomatic/mc3.v0.2.8.PUBLIC.BRCA.func.rmdup.maf",sep="\t")
#filter high impact somatic
somaticSNVImpact<-somaticSNV[-intersect(grep("tolerated",somaticSNV$SIFT), grep("benign",somaticSNV$PolyPhen)),]

CEUsample<-read.table("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/PCACombindGermlineSomaticMutation/Europen_IndividualID_in_TCGA_BRCA_Final.txt",h=F)
clinical<-read.csv("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/PCACombindGermlineSomaticMutation/PanCan_ClinicalData_V4_20170428_BRCA.txt",h=T,sep="\t")
CEUsomaticSNV<-somaticSNV[which(somaticSNV$Tumor_Sample_Barcode%in%CEUsample[,1]),]
CEUgermSNV<-germSNV[which(germSNV$Tumor_Sample_Barcode%in%CEUsample[,1]),]

#pathways 
load("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/Pathway/Pathway.RData")
names(hallmarkpathway)<-gsub("HALLMARK_","",names(hallmarkpathway))
BiocartPathway<-BiocartPathway[-which(lapply(BiocartPathway,length)<5)]
CorePathway<-CorePathway[-which(lapply(CorePathway,length)<5)]
#####################################Pathway Interaction Analysis#########################################

pathInteract<-function(germSNV,somaticSNV,samplename,pathLabel,pathname,pathfile,pvalue = c(0.05, 0.01), fontSize = 1){
  #construct pathway mutation status matrix
  GermlinePathaltMat<-matrix(data="",ncol=length(samplename),nrow=length(pathname),dimnames=list(pathname,samplename))
  SomaticPathaltMat<-matrix(data="",ncol=length(samplename),nrow=length(pathname),dimnames=list(pathname,samplename))
  
  for(i in 1:length(pathname)){
    pathgene<-pathfile[[i]]
    #germline
    if(length(which(germSNV$Hugo_Symbol%in%pathgene))!=0){
      germlinepathgenemat<-germSNV[which(germSNV$Hugo_Symbol%in%pathgene),]
      GermlineMutSpl<-intersect(germlinepathgenemat$Tumor_Sample_Barcode,CEUsample[,1])
      if(length(GermlineMutSpl)!=0){
        GermlinePathaltMat[pathname[i],GermlineMutSpl]<-"GermlineMutations"
      }
    }
    #somatic
    if(length(which(somaticSNV$Hugo_Symbol%in%pathgene))!=0){
      somaticpathgenemat<-somaticSNV[which(somaticSNV$Hugo_Symbol%in%pathgene),]
      SomaticMutSpl<-intersect(somaticpathgenemat$Tumor_Sample_Barcode,CEUsample[,1])
      SomaticPathaltMat[pathname[i],SomaticMutSpl]<-"SomaticMutations"
    }
    print(i)
  }
  
  
  write.csv(GermlinePathaltMat,paste0(pathLabel,"_GermlinePathaltMat.csv"))
  write.csv(SomaticPathaltMat,paste0(pathLabel,"SomaticPathaltMat.csv"))
  #oncoprinter: draw pathway mutation map
  require(ComplexHeatmap)
  col = c("GermlineMutations" = "forestgreen", "SomaticMutations" = "slateblue4")
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    GermlineMutations = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.65, "mm"), h-unit(1.65, "mm"), gp = gpar(fill = "forestgreen", col = NA))
    },
    SomaticMutations = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.65, "mm"), h-unit(1.65, "mm"), gp = gpar(fill = "slateblue4", col = NA))
    }
  )
  
  #highimpact<-order(apply(GermlinePathaltMat,1,function(x)length(grep("GermlineMutations",x))),decreasing=T)[1:50]
  #GermlinePathaltMat<-GermlinePathaltMat[highimpact,]
  #highimpact<-order(apply(SomaticPathaltMat,1,function(x)length(grep("SomaticMutations",x))),decreasing=T)[1:50]
  #SomaticPathaltMat<-SomaticPathaltMat[highimpact,]
  
  
  png(paste0(pathLabel,"_PCA_pathVargermline_Oncoprinter.png"),height=2000,width=2000,res=96)
  print(oncoPrint(GermlinePathaltMat,  alter_fun = alter_fun,col = col,get_type = function(x) strsplit(x, ";")[[1]], column_title =  paste0("Breast Cancer Germline Pathogenic Mutations (",pathLabel,")"),
                  heatmap_legend_param = list(title = "Alternations", at = c("GermlineMutations" , "SomaticMutations"), 
                                              labels = c("GermlineMutations" , "SomaticMutations")
                  )))
  dev.off()
  
  png(paste0(pathLabel,"_PCA_pathVarSomatic_Oncoprinter.png"),height=1000,width=2000,res=96)
  print(oncoPrint(SomaticPathaltMat,  alter_fun = alter_fun,col = col,get_type = function(x) strsplit(x, ";")[[1]], column_title = paste0("Breast Cancer Somatic Mutations (",pathLabel,")"),
                  heatmap_legend_param = list(title = "Alternations", at = c("GermlineMutations" , "SomaticMutations"),
                                              labels = c("GermlineMutations" , "SomaticMutations")
                  )))
  dev.off()
  
  ##interaction analysis
  
  GermlinePathaltMat1<-GermlinePathaltMat
  GermlinePathaltMat1[which(GermlinePathaltMat1=="")]=0
  GermlinePathaltMat1[which(GermlinePathaltMat1=="GermlineMutations")]=1
  
  SomaticPathaltMat1<-SomaticPathaltMat
  SomaticPathaltMat1[which(SomaticPathaltMat1=="")]=0
  SomaticPathaltMat1[which(SomaticPathaltMat1=="SomaticMutations")]=1
  
  interactions = sapply(pathname, function(i) sapply(pathname, function(j) {f<- try(fisher.test(GermlinePathaltMat1[i,], SomaticPathaltMat1[j,]), silent=TRUE); if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.value),log10(f$p.value))} ))
  oddsRatio <- oddsGenes <- sapply(pathname, function(i) sapply(pathname, function(j) {f<- try(fisher.test(GermlinePathaltMat1[i,], SomaticPathaltMat1[j,]), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = pathname
  
  #interaction figure
  #Source code borrowed from: https://www.nature.com/articles/ncomms6901
  interactions[10^-abs(interactions) > max(pvalue)] = 0
  diag(interactions) <- 0
  m <- nrow(interactions)
  n <- ncol(interactions)
  interactions[interactions < -4] = -4
  interactions[interactions > 4] = 4
  interactions[lower.tri(x = interactions)] = NA
  
  png(paste0(pathLabel,"_somatic_germline_interaction.png"),height=1000,width=1000)  
  par(bty="n", mgp = c(2,.5,0), mar = c(2, 20, 20, 5)+.1, las=2, tcl=-.33)
  image(x=1:n, y=1:m, interactions, col=RColorBrewer::brewer.pal(9,"PiYG"),
        breaks = c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n",
        xlab="Germline",ylab="Somatic", xlim=c(0, n+4), ylim=c(0, n+4))
  abline(h=0:n+.5, col="white", lwd=.5)
  abline(v=0:n+.5, col="white", lwd=.5)
  
  mtext(side = 2, at = 1:m, text = colnames(interactions), cex = fontSize, font = 2)
  mtext(side = 3, at = 1:n, text = colnames(interactions), las = 2, line = -2, cex = fontSize, font = 2)
  
  w = arrayInd(which(10^-abs(interactions) < min(pvalue)), rep(m,2))
  points(w, pch="*", col="black")
  w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
  points(w, pch=".", col="black")
  image(y = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 8), x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col = RColorBrewer::brewer.pal(8,"PiYG"), add=TRUE)
  atLims = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 7)
  axis(side = 4, at = atLims,  tcl=-.15, labels =c(3:1, 0, 1:3), las=1, lwd=.5)
  mtext(side=4, at = median(atLims), "log10 (p-value)", las=3, cex = 0.9, line = 3, font = 2)
  
  par(xpd=NA)
  text(x=n+2.2, y= max(atLims)+1.2, "Co-occurance", pos=4, cex = 0.9, font = 2)
  text(x=n+2.2, y = min(atLims)-1.2, "Exclusive", pos=4, cex = 0.9, font = 2)
  points(x = n+1, y = 0.2*n, pch = "*", cex = 2)
  text(x = n+1, y = 0.2*n, paste0(" p < ", min(pvalue)), pos=4, cex = 0.9, font = 2)
  points(x = n+1, y = 0.1*n, pch = ".", cex = 2)
  text(x = n+1, y = 0.1*n, paste0("p < ", max(pvalue)), pos=4, cex = 0.9)
  dev.off()
}



pathInteract(germSNV=germSNV,somaticSNV=somaticSNVImpact,samplename=CEUsample[,1],pathLabel="HallMark",pathname=names(hallmarkpathway),pathfile=hallmarkpathway)

pathInteract(germSNV=germSNV,somaticSNV=somaticSNVImpact,samplename=CEUsample[,1],pathLabel="Biocart",pathname=names(BiocartPathway),pathfile=BiocartPathway)

pathInteract(germSNV=germSNV,somaticSNV=somaticSNVImpact,samplename=CEUsample[,1],pathLabel="CorePathways",pathname=names(CorePathway),pathfile=CorePathway)










