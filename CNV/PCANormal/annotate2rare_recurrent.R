#this script annotate PCA normal CNV to the rare/recurrent DGV database

##################################functions#####################################################

geneRanges <- function(db, column="ENTREZID"){
  #genes
  g0 <- genes(db, columns=column)
  g<-g0[which(seqnames(g0)%in%paste0("chr",c(1:22,"X","Y")))]
  col <- mcols(g)[[column]]
  genes <- g[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  
  #CDS
  c0 <- cds(db, columns=column)
  c<-c0[which(seqnames(c0)%in%paste0("chr",c(1:22,"X","Y")))]
  cdscol <- mcols(c)[[column]]
  cds <- granges(c)[rep(seq_along(c), elementNROWS(cdscol))]
  mcols(cds)[[column]] <- as.character(unlist(cdscol))
  
  #exon
  e0 <- exons(db, columns=column)
  e<-e0[which(seqnames(e0)%in%paste0("chr",c(1:22,"X","Y")))]
  exoncol <- mcols(e)[[column]]
  exon <- granges(e)[rep(seq_along(e), elementNROWS(exoncol))]
  mcols(exon)[[column]] <- as.character(unlist(exoncol))
  
  #promoters
  p0 <- promoters(db, upstream=2000, downstream=200,columns=c("gene_id","tx_name"))
  p<-p0[which(seqnames(p0)%in%paste0("chr",c(1:22,"X","Y")))]
  procol <- mcols(p)[["gene_id"]]
  proms <- granges(p)[rep(seq_along(p), elementNROWS(procol))]
  mcols(proms)[["gene_id"]] <- as.character(unlist(procol))
  
  anno=list(genes=genes,cds=cds,exons=exon,promoters=proms)
  return(anno)
}


CNVFilter <- function(query,subject,overlap=0.7,label,column="ENTREZID", ...){
    #Compute percent overlap with DGV and filter the hits:
    hits <- findOverlaps(query, subject)
	  overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
	  percentOverlap <- width(overlaps) / width(subject[subjectHits(hits)])
	  olaps <- hits[percentOverlap > overlap]
    f1 <- factor(subjectHits(olaps),levels=seq_len(subjectLength(olaps)))
    cnvMatRecurrent<-as.data.frame(subject)[unique(f1),]
    cnvMatRare<-as.data.frame(subject)[-as.numeric((unique(f1))),]
    fwrite(as.data.frame(cnvMatRecurrent),paste0(label,"recurrent_overlaped_DGV.txt"),quote=F,row.names = F,sep="\t")
    fwrite(as.data.frame(cnvMatRare),paste0(label,"rare_overlaped_DGV.txt"),quote=F,row.names = F,sep="\t")
    return(list(cnvMatRecurrent=cnvMatRecurrent,cnvMatRare=cnvMatRare))
}

CNVAnno <- function(query,subject,overlap=0.1,label,type,column="ENTREZID", ...){
    hits <- findOverlaps(query, subject)
    overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(subject[subjectHits(hits)])
    olaps <- hits[percentOverlap > overlap]
    f1 <- factor(subjectHits(olaps),levels=seq_len(subjectLength(olaps)))
    cnvMat<-as.data.frame(subject) 
    
    hitedgenelist<-splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
    hitedgenelistSep<-unlist(lapply(hitedgenelist,function(x)paste(unique(x),collapse = "|")))
    cnvMatanno<-cbind(cnvMat,hitedgenelistSep)
    lens=dim(cnvMatanno)[2]
    
    require(foreach)
    require(doParallel)
    registerDoParallel(cores=30)
    
    cnvGeneMat=foreach(i=c(1:dim(cnvMatanno)[1]),.combine = "rbind")%dopar%{
      print(i)
      iterm<-as.character(as.matrix(cnvMatanno[i,]))
      genes=unlist(strsplit(iterm[lens],split="\\|"))
      cbind(matrix(rep(as.character(iterm[1:lens]),length(genes)),ncol=lens,byrow=T),genes)
    }   

    cnvGeneMat<-as.data.frame(cnvGeneMat[,-5])
    #colnames(cnvGeneMat)<-c("Chr","Start","End","VariantType","SVStatus","NumberofVariants","NumberofStudies","NumberofSamples","TotalSamples","Frequency","African","Asia","European","Mexican","MiddleEast","NativeAmerican","NorthAmerican","Oceania","SourthAmerican","Turkish","Admixed","Unkown","Info","Gene")
    colnames(cnvGeneMat)<-c("Chr","start","end","length","sample","NumbProbe","log2Ratio","SampleType","CNVStatus","Anno","Genes") 
    cnvGeneMat<-as.data.frame(cnvGeneMat)
    if(any(cnvGeneMat$Genes=="NA")){
      cnvGeneMat<-cnvGeneMat[-which(cnvGeneMat$Genes=="NA"),]
    }
    cnvGeneMat<- apply(cnvGeneMat,2,as.character)
    
    if(type=="genes"||type=="cds"||type=="exon"){
      cnvGeneMat<-cbind(cnvGeneMat,lapply(mget(x=cnvGeneMat[,(lens)],envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
      colnames(cnvGeneMat)[lens+1]<-"gene_id"
    }
      
    if(type=="promoters"){
      require(annotate)
      cnvGeneMat<-cbind(cnvGeneMat[,1:(lens-1)],unlist(lapply(getSYMBOL(cnvGeneMat[,lens], data='org.Hs.eg'),function(x)paste(unique(x),collapse = "|"))),cnvGeneMat[,(lens)])
      colnames(cnvGeneMat)[c((lens),(lens+1))]<-c("Gene","gene_id")
    }
    
    fwrite(cnvMatanno,paste0(label,"cnvMatanno_",type,"level.txt"),sep="\t",row.names=F)
    fwrite(as.data.frame(cnvGeneMat),paste0(label,"cnvGeneMat_",type,"level.txt"),quote=F,sep="\t",row.names=F)
    return(list(cnvMatanno=cnvMatanno,cnvGeneMat=cnvGeneMat,lens=lens))
}

#parameters
#args = commandArgs(trailingOnly=TRUE)
#Chr	Start	End	VariantType	SVStatus	NumberofVariants	NumberofStudies	NumberofSamples	TotalSamples	Frequency	African	Asia	European	Mexican	MiddleEast	NativeAmerican	NorthAmerican	Oceania	SourthAmerican	Turkish	Admixed	Unkown	Info

cnvgainfile="PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_gain.txt"#args[1]
cnvlossfile="PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_loss.txt"
#DGV file
dgvRecurrentGainFile<-"/home/tao/taoqing/references/DGV/DGV.GoldStanderd.matrix.March2016.50percent_HC_recurrent_gain.txt"
dgvRecurrentLossFile<-"/home/tao/taoqing/references/DGV/DGV.GoldStanderd.matrix.March2016.50percent_HC_recurrent_loss.txt"
dgvRareGainFile<-"/home/tao/taoqing/references/DGV/DGV.GoldStanderd.matrix.March2016.50percent_HC_rare_gain.txt"
dgvRareLossFile<-"/home/tao/taoqing/references/DGV/DGV.GoldStanderd.matrix.March2016.50percent_HC_rare_loss.txt"

library(data.table)
cnvGainmat<-as.data.frame(fread(cnvgainfile))
cnvLossmat<-as.data.frame(fread(cnvlossfile))
dgvRecurrentGainMat<-as.data.frame(fread(dgvRecurrentGainFile))
dgvRecurrentLossMat<-as.data.frame(fread(dgvRecurrentLossFile))
dgvRareGainMat<-as.data.frame(fread(dgvRareGainFile))
dgvRareLossMat<-as.data.frame(fread(dgvRareLossFile))

###################################CNA annotation##################################################
#prepare range data
require(GenomicRanges)
cnvgain = makeGRangesFromDataFrame(cnvGainmat,keep.extra.columns=T)
cnvloss = makeGRangesFromDataFrame(cnvLossmat,keep.extra.columns=T)
dgvRecurrentGain<-makeGRangesFromDataFrame(dgvRecurrentGainMat,keep.extra.columns=T)
dgvRecurrentLoss<-makeGRangesFromDataFrame(dgvRecurrentLossMat,keep.extra.columns=T)
dgvRareGain<-makeGRangesFromDataFrame(dgvRareGainMat,keep.extra.columns=T)
dgvRareLoss<-makeGRangesFromDataFrame(dgvRareLossMat,keep.extra.columns=T)

require(Homo.sapiens)
annotation = geneRanges(Homo.sapiens, column="SYMBOL")#gene annotation

#########################################Recurrent Gain########################################
#genelevel
gns=annotation$genes
PCAGDV_Gain_Overlap = CNVFilter(dgvRecurrentGain,cnvgain,label="PCAnormal_cnv_gain_",column="SYMBOL",overlap=0.7)#query=dgvRecurrentGain;subject=cnvgain;geneanno=gns;column="SYMBOL";overlap=0.7
PCAGDV_RecGain_OverlapMat<-makeGRangesFromDataFrame(PCAGDV_Gain_Overlap$cnvMatRecurrent,keep.extra.columns=T)
PCAGDV_RareGain_OverlapMat<-makeGRangesFromDataFrame(PCAGDV_Gain_Overlap$cnvMatRare,keep.extra.columns=T)

PCAGDV_Loss_Overlap = CNVFilter(dgvRecurrentLoss,cnvloss,label="PCAnormal_cnv_loss_",column="SYMBOL",overlap=0.7)#query=dgvRecurrentGain;subject=cnvgain;geneanno=gns;column="SYMBOL";overlap=0.7
PCAGDV_RecLoss_OverlapMat<-makeGRangesFromDataFrame(PCAGDV_Loss_Overlap$cnvMatRecurrent,keep.extra.columns=T)
PCAGDV_RareLoss_OverlapMat<-makeGRangesFromDataFrame(PCAGDV_Loss_Overlap$cnvMatRare,keep.extra.columns=T)



PCAGDV_RecGain_OverlapMatGeneAnno <- CNVAnno(gns,PCAGDV_RecGain_OverlapMat,label="PCANormal_CNV_Recurrent_Gain_",overlap=0.1,column="SYMBOL",type="genes")#query=gns;subject=PCAGDV_RecGain_OverlapMat;overlap=0.1;column="SYMBOL"
PCAGDV_RecLoss_OverlapMatGeneAnno <- CNVAnno(gns,PCAGDV_RecLoss_OverlapMat,label="PCANormal_CNV_Recurrent_Loss_",overlap=0.1,column="SYMBOL",type="genes")
PCAGDV_RareGain_OverlapMatGeneAnno <- CNVAnno(gns,PCAGDV_RareGain_OverlapMat,label="PCANormal_CNV_Rare_Gain_",overlap=0.1,column="SYMBOL",type="genes")
PCAGDV_RareLoss_OverlapMatGeneAnno <- CNVAnno(gns,PCAGDV_RareLoss_OverlapMat,label="PCANormal_CNV_Rare_Loss_",overlap=0.1,column="SYMBOL",type="genes")


PCAGDV_RecGain_OverlapMatCDSAnno <- CNVAnno(annotation$cds,PCAGDV_RecGain_OverlapMat,label="PCANormal_CNV_Recurrent_Gain_",overlap=0.1,column="SYMBOL",type="cds")
PCAGDV_RecLoss_OverlapMatCDSAnno <- CNVAnno(annotation$cds,PCAGDV_RecLoss_OverlapMat,label="PCANormal_CNV_Recurrent_Loss_",overlap=0.1,column="SYMBOL",type="cds")
PCAGDV_RareGain_OverlapMatCDSAnno <- CNVAnno(annotation$cds,PCAGDV_RareGain_OverlapMat,label="PCANormal_CNV_Rare_Gain_",overlap=0.1,column="SYMBOL",type="cds")
PCAGDV_RareLoss_OverlapMatCDSAnno <- CNVAnno(annotation$cds,PCAGDV_RareLoss_OverlapMat,label="PCANormal_CNV_Rare_Loss_",overlap=0.1,column="SYMBOL",type="cds")

PCAGDV_RecGain_OverlapMatExonAnno <- CNVAnno(annotation$exons,PCAGDV_RecGain_OverlapMat,label="PCANormal_CNV_Recurrent_Gain_",overlap=0.1,column="SYMBOL",type="exon")
PCAGDV_RecLoss_OverlapMatExonAnno <- CNVAnno(annotation$exons,PCAGDV_RecLoss_OverlapMat,label="PCANormal_CNV_Recurrent_Loss_",overlap=0.1,column="SYMBOL",type="exon")
PCAGDV_RareGain_OverlapMatExonAnno <- CNVAnno(annotation$exons,PCAGDV_RareGain_OverlapMat,label="PCANormal_CNV_Rare_Gain_",overlap=0.1,column="SYMBOL",type="exon")
PCAGDV_RareLoss_OverlapMatExonAnno <- CNVAnno(annotation$exons,PCAGDV_RareLoss_OverlapMat,label="PCANormal_CNV_Rare_Loss_",overlap=0.1,column="SYMBOL",type="exon")

PCAGDV_RecGain_OverlapMatPromotersAnno <- CNVAnno(annotation$promoters,PCAGDV_RecGain_OverlapMat,label="PCANormal_CNV_Recurrent_Gain_",overlap=0.1,column="gene_id",type="promoters")
PCAGDV_RecLoss_OverlapMatPromotersAnno <- CNVAnno(annotation$promoters,PCAGDV_RecLoss_OverlapMat,label="PCANormal_CNV_Recurrent_Loss_",overlap=0.1,column="gene_id",type="promoters")
PCAGDV_RareGain_OverlapMatPromotersAnno <- CNVAnno(annotation$promoters,PCAGDV_RareGain_OverlapMat,label="PCANormal_CNV_Rare_Gain_",overlap=0.1,column="gene_id",type="promoters")
PCAGDV_RareLoss_OverlapMatPromotersAnno <- CNVAnno(annotation$promoters,PCAGDV_RareLoss_OverlapMat,label="PCANormal_CNV_Rare_Loss_",overlap=0.1,column="gene_id",type="promoters")
#query=annotation$promoters;subject=PCAGDV_RareLoss_OverlapMat;label="PCANormal_CNV_Rare_Loss_";overlap=0.1;column="SYMBOL";type="promoters"
#"SYMBOL","gene_id",annotation$cds,annotation$exons,annotation$promoters
save.image("PCANormal_CNV_Classification.RData")


