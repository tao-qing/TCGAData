#this script annotated the DGV database into four genomic level: gene, exon, CDS, promoters, with a specified overlap 

##################################functions#####################################################
geneRangesold <- function(db, column="ENTREZID"){
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}


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

CNVanno <- function(query, subject,overlap=0.10, column="ENTREZID", ...){
  hits <- findOverlaps(query, subject)
  #Compute percent overlap and filter the hits:
  overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(subject[subjectHits(hits)])
  olaps <- hits[percentOverlap > overlap]
  
  f1 <- factor(subjectHits(olaps),levels=seq_len(subjectLength(olaps)))
  cnvMat<-as.data.frame(subject)
  
  hitedgenelist<-splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  hitedgenelistSep<-unlist(lapply(hitedgenelist,function(x)paste(unique(x),collapse = "|")))
  cnvMatanno<-cbind(cnvMat,hitedgenelistSep)
  if(any(hitedgenelistSep=="")){
    cnvMatanno<-cnvMatanno[-which(hitedgenelistSep==""),]
  }
  lens=dim(cnvMatanno)[2]
  
  #for(i in 1:dim(cnvMatanno)[1]){
  #  iterm<-as.character(as.matrix(cnvMatanno[i,]))
  #  genes=unlist(strsplit(iterm[lens],split="\\|"))
  #  if(i==1){
  #    cnvGeneMat<-cbind(matrix(rep(as.character(iterm[1:lens]),length(genes)),ncol=lens,byrow=T),genes)
  #    # mergegene<-genes
  #  }else{
  #    cnvGeneMat<-rbind(cnvGeneMat,cbind(matrix(rep(as.character(iterm[1:lens]),length(genes)),ncol=lens,byrow=T),genes))
  #    #  mergegene<-c(mergegene,genes)
  # }
  #  print (i)
  #}
  
  require(foreach)
  require(doParallel)
  registerDoParallel(cores=30)
  
  cnvGeneMat=foreach(i=c(1:dim(cnvMatanno)[1]),.combine = "rbind")%dopar%{
    print(i)
    iterm<-as.character(as.matrix(cnvMatanno[i,]))
    genes=unlist(strsplit(iterm[lens],split="\\|"))
    cbind(matrix(rep(as.character(iterm[1:lens]),length(genes)),ncol=lens,byrow=T),genes)
  }
  
  
  #cnvGeneMat[,9]<-mergegene
  cnvGeneMat<-as.data.frame(cnvGeneMat[,-5])
  #colnames(cnvGeneMat)<-c("Chr","Start","End","VariantType","SVStatus","NumberofVariants","NumberofStudies","NumberofSamples","TotalSamples","Frequency","African","Asia","European","Mexican","MiddleEast","NativeAmerican","NorthAmerican","Oceania","SourthAmerican","Turkish","Admixed","Unkown","Info","Gene")
  colnames(cnvGeneMat)<-c("Chr","start","end","length","sample","NumbProbe","log2Ratio","SampleType","CNVStatus","Anno","Genes") 
  if(any(cnvGeneMat$Genes=="NA")){
    cnvGeneMat<-cnvGeneMat[-which(cnvGeneMat$Genes=="NA"),]
  }
  cnvGeneMat<- apply(cnvGeneMat,2,as.character)
  return(list(cnvMatanno=cnvMatanno,cnvGeneMat=cnvGeneMat,lens=lens))
}



cnvfile="PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal.txt"#args[1]

library(data.table)
rawmat<-as.data.frame(fread(cnvfile))
###################################CNA annotation##################################################
require(GenomicRanges)
cnv = makeGRangesFromDataFrame(rawmat,keep.extra.columns=T)
#cnv = GRanges(c("chr1", "chr2"), IRanges(94312388,244006886))

require(Homo.sapiens)
annotation = geneRanges(Homo.sapiens, column="SYMBOL")

#genelevel
gns=annotation$genes
geneOut = CNVanno(gns, cnv, "SYMBOL",overlap=0.1)

#cds level
cds=annotation$cds
cdsOut = CNVanno(cds, cnv, "SYMBOL",overlap=0.1)

#exon level
exon=annotation$exons
exonsOut = CNVanno(exon, cnv, "SYMBOL",overlap=0.1)

#promoters level
proms=annotation$promoters
promsOut = CNVanno(proms, cnv, "gene_id",overlap=0.1)

########add gene id######################
gns_cnvMatanno<-geneOut$cnvMatanno
gns_cnvGeneMat<-geneOut$cnvGeneMat
lens<-geneOut$lens
gns_cnvGeneMat<-cbind(gns_cnvGeneMat,lapply(mget(as.character(gns_cnvGeneMat[,(lens)]),envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
colnames(gns_cnvGeneMat)[lens+1]<-"gene_id"


cds_cnvMatanno<-cdsOut$cnvMatanno
cds_cnvGeneMat<-cdsOut$cnvGeneMat
lens<-cdsOut$lens
cds_cnvGeneMat<-cbind(cds_cnvGeneMat,lapply(mget(as.character(cds_cnvGeneMat[,(lens)]),envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
colnames(cds_cnvGeneMat)[lens+1]<-"gene_id"


exon_cnvMatanno<-exonsOut$cnvMatanno
exon_cnvGeneMat<-exonsOut$cnvGeneMat
lens<-exonsOut$lens
exon_cnvGeneMat<-cbind(exon_cnvGeneMat,lapply(mget(as.character(exon_cnvGeneMat[,(lens)]),envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
colnames(exon_cnvGeneMat)[lens+1]<-"gene_id"


require(annotate)
proms_cnvMatanno<-promsOut$cnvMatanno
proms_cnvGeneMat<-promsOut$cnvGeneMat
lens<-promsOut$lens
proms_cnvGeneMat<-cbind(proms_cnvGeneMat[,1:lens-1],lapply(getSYMBOL(proms_cnvGeneMat[,(lens)], data='org.Hs.eg'),function(x)paste(unique(x),collapse = "|")),proms_cnvGeneMat[,(lens)])
colnames(proms_cnvGeneMat)[c((lens),(lens+1))]<-c("Gene","gene_id")


fwrite(gns_cnvMatanno,"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvMatanno_genelevel.txt",quote=F,sep="\t",row.names=F)
fwrite(as.data.frame(gns_cnvGeneMat),"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvGeneMat_genelecel.txt",quote=F,sep="\t",row.names=F)

fwrite(cds_cnvMatanno,"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvMatanno_cdslevel.txt",quote=F,sep="\t",row.names=F)
fwrite(as.data.frame(cds_cnvGeneMat),"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvGeneMat_cdslevel.txt",quote=F,sep="\t",row.names=F)

fwrite(exon_cnvMatanno,"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvMatanno_exonlevel.txt",quote=F,sep="\t",row.names=F)
fwrite(as.data.frame(exon_cnvGeneMat),"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvGeneMat_exonlevel.txt",quote=F,sep="\t",row.names=F)

fwrite(proms_cnvMatanno,"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvMatanno_promoterlevel.txt",quote=F,sep="\t",row.names=F)
fwrite(as.data.frame(proms_cnvGeneMat),"PCA_germline_cnv_filtered_segmean0.1_probe10_nosex_normal_cnvGeneMat_promoterlevel.txt",quote=F,sep="\t",row.names=F)

save.image("DGV.GoldStanderd.Annotation.RData")
print (paste0("Copy number annotation of  Finshed!"))

