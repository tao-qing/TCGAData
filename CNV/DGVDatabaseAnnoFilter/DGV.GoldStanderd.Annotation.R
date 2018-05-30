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
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))

    #CDS
    c <- cds(db, columns=column)
    cdscol <- mcols(c)[[column]]
    cds <- granges(c)[rep(seq_along(c), elementNROWS(cdscol))]
    mcols(cds)[[column]] <- as.character(unlist(cdscol))

    #exon
    e <- exons(db, columns=column)
    exoncol <- mcols(e)[[column]]
    exon <- granges(e)[rep(seq_along(e), elementNROWS(exoncol))]
    mcols(exon)[[column]] <- as.character(unlist(exoncol))

    #promoters
    p <- promoters(db, upstream=2000, downstream=200,columns=c("gene_id","tx_name"))
    procol <- mcols(p)[["gene_id"]]
    proms <- granges(p)[rep(seq_along(p), elementNROWS(procol))]
    mcols(proms)[["gene_id"]] <- as.character(unlist(procol))

    anno=list(genes=genes,cds=cds,exons=exon,promoters=proms)
    return(anno)
  }



CNVanno <- function(query, subject,overlap=0.5, column="ENTREZID", ...){
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
	lens=dim(cnvMatanno)[2]

    for(i in 1:dim(cnvMatanno)[1]){
      iterm<-as.character(as.matrix(cnvMatanno[i,]))
      genes=unlist(strsplit(iterm[lens],split="\\|"))
      if(i==1){
        cnvGeneMat<-cbind(matrix(rep(as.character(iterm[1:lens]),length(genes)),ncol=lens,byrow=T),genes)
       # mergegene<-genes
      }else{
        cnvGeneMat<-rbind(cnvGeneMat,cbind(matrix(rep(as.character(iterm[1:lens]),length(genes)),ncol=lens,byrow=T),genes))
      #  mergegene<-c(mergegene,genes)
      }
	print (i)
    }
     
    #cnvGeneMat[,9]<-mergegene
    cnvGeneMat<-as.data.frame(cnvGeneMat[,-5])
    colnames(cnvGeneMat)<-c("Chr","Start","End","VariantType","SVStatus","NumberofVariants","NumberofStudies","NumberofSamples","TotalSamples","Frequency","African","Asia","European","Mexican","MiddleEast","NativeAmerican","NorthAmerican","Oceania","SourthAmerican","Turkish","Admixed","Unkown","Info","Gene")
    if(any(cnvGeneMat$Genes=="NA")){
      cnvGeneMat<-cnvGeneMat[-which(cnvGeneMat$Genes=="NA"),]
    }
    cnvGeneMat<- apply(cnvGeneMat,2,as.character)
    return(list(cnvMatanno=cnvMatanno,cnvGeneMat=cnvGeneMat))
}

#parameters
#args = commandArgs(trailingOnly=TRUE)
#Chr	Start	End	VariantType	SVStatus	NumberofVariants	NumberofStudies	NumberofSamples	TotalSamples	Frequency	African	Asia	European	Mexican	MiddleEast	NativeAmerican	NorthAmerican	Oceania	SourthAmerican	Turkish	Admixed	Unkown	Info

cnvfile="DGV.GoldStander.matrix.March2016.50percent.txt"#args[1]

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
geneOut = CNVanno(gns, cnv, "SYMBOL",overlap=0.4)

#cds level
cds=annotation$cds
cdsOut = CNVanno(cds, cnv, "SYMBOL",overlap=1)

#exon level
exons=annotation$exons
exonsOut = CNVanno(exons, cnv, "SYMBOL",overlap=1)

#promoters level
proms=annotation$promoters
promsOut = CNVanno(proms, cnv, "gene_id",overlap=0.4)



#rawmat[,2]<-paste0("chr",rawmat[,2])
gns_cnvMatanno<-geneOut$cnvMatanno
gns_cnvGeneMat<-geneOut$cnvGeneMat
gns_cnvGeneMat<-cbind(gns_cnvGeneMat,lapply(mget(x=gns_cnvGeneMat[,(lens+1)],envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
colnames(gns_cnvGeneMat)[lens+2]<-"gene_id"


cds_cnvMatanno<-cdsOut$cnvMatanno
cds_cnvGeneMat<-cdsOut$cnvGeneMat
cds_cnvGeneMat<-cbind(cds_cnvGeneMat,lapply(mget(x=cds_cnvGeneMat[,(lens+1)],envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
colnames(cds_cnvGeneMat)[lens+2]<-"gene_id"


exon_cnvMatanno<-exonsOut$cnvMatanno
exon_cnvGeneMat<-exonsOut$cnvGeneMat
exon_cnvGeneMat<-cbind(exon_cnvGeneMat,lapply(mget(x=exon_cnvGeneMat[,(lens+1)],envir=org.Hs.egALIAS2EG),function(x)paste(unique(x),collapse = "|")))
colnames(exon_cnvGeneMat)[lens+2]<-"gene_id"


require(annotate)
proms_cnvMatanno<-promsOut$cnvMatanno
proms_cnvGeneMat<-promsOut$cnvGeneMat
proms_cnvGeneMat<-cbind(proms_cnvGeneMat[,1:lens],lapply(getSYMBOL(proms_cnvGeneMat[,1:(lens+1)], data='org.Hs.eg'),function(x)paste(unique(x),collapse = "|")),proms_cnvGeneMat[,(lens+1)])
colnames(proms_cnvGeneMat)[c((lens+1),(lens+2)]<-c("Gene","gene_id")


write.table(gns_cnvMatanno,"DGV_GoldStanderd_cnvMatanno_genelevel.txt",sep="\t",row.names=F)
write.table(gns_cnvGeneMat,"DGV_GoldStanderd_cnvGeneMat_genelecel.txt",quote=F,sep="\t",row.names=F)

write.table(cds_cnvMatanno,"DGV_GoldStanderd_cnvMatanno_cdslevel.txt",sep="\t",row.names=F)
write.table(cds_cnvGeneMat,"DGV_GoldStanderd_cnvGeneMat_cdslevel.txt",quote=F,sep="\t",row.names=F)

write.table(exon_cnvMatanno,"DGV_GoldStanderd_cnvMatanno_exonlevel.txt",sep="\t",row.names=F)
write.table(exon_cnvGeneMat,"DGV_GoldStanderd_cnvGeneMat_exonlevel.txt",quote=F,sep="\t",row.names=F)

write.table(proms_cnvMatanno,"DGV_GoldStanderd_cnvMatanno_promoterlevel.txt",sep="\t",row.names=F)
write.table(proms_cnvGeneMat,"DGV_GoldStanderd_cnvGeneMat_promoterlevel.txt",quote=F,sep="\t",row.names=F)


print (paste0("Copy number annotation of  Finshed!"))

