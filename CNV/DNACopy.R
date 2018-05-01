###################################functions#####################################################
geneRanges <- function(db, column="ENTREZID"){
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }

CNVanno <- function(query, subject, column="ENTREZID", ...){
    olaps <- findOverlaps(query, subject)
    f1 <- factor(subjectHits(olaps),levels=seq_len(subjectLength(olaps)))
    cnvMat<-as.data.frame(subject)
    
    hitedgenelist<-splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
    hitedgenelistSep<-unlist(lapply(hitedgenelist,function(x)paste(x,collapse = "|")))
    cnvMatanno<-cbind(cnvMat,hitedgenelistSep)

    for(i in 1:dim(cnvMatanno)[1]){
      iterm<-as.character(as.matrix(cnvMatanno[i,]))
      genes=unlist(strsplit(iterm[9],split="\\|"))
      if(i==1){
        cnvGeneMat<-matrix(rep(iterm,length(genes)),ncol=9,byrow=T)
        mergegene<-genes
      }else{
        cnvGeneMat<-rbind(cnvGeneMat,matrix(rep(cnvMatanno[i,],length(genes)),ncol=9,byrow=T))
        mergegene<-c(mergegene,genes)
      }
    }
    cnvGeneMat[,9]<-mergegene
    cnvGeneMat<-as.data.frame(cnvGeneMat[,-5])
    colnames(cnvGeneMat)<-c("chr","start","end","length","sample","NumbProbe","log2Ratio","Genes")
    if(any(cnvGeneMat$Genes=="NA")){
      cnvGeneMat<-cnvGeneMat[-which(cnvGeneMat$Genes=="NA"),]
    }
    cnvGeneMat<- apply(cnvGeneMat,2,as.character)
    return(list(cnvMatanno=cnvMatanno,cnvGeneMat=cnvGeneMat))
}

#parameters
args = commandArgs(trailingOnly=TRUE)

cnvfile=args[1]
spl_name=args[2]

library(data.table)
rawmat<-as.data.frame(fread(cnvfile))
matrmna<-rawmat[-which(is.na(rawmat[,6])),]
matrmna[,6]<-round(log2((matrmna[,6])/2),digits=4)
###################################Sigmentation analysis###########################################
require(DNAcopy)
CNA.object <- CNA(as.numeric(as.matrix(matrmna[,6])),matrmna$CHROM,as.numeric(as.matrix(matrmna$POS)),data.type="logratio",sampleid=spl_name)
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

png(paste0(spl_name,"_cnv.png"),width=2000,height=500)
plot(segment.smoothed.CNA.object, plot.type="s",xmaploc=T)
dev.off()

###################################CNA annotation##################################################
require(GenomicRanges)
df = segment.smoothed.CNA.object$output
cnv = makeGRangesFromDataFrame(df,keep.extra.columns=T)
#cnv = GRanges(c("chr1", "chr2"), IRanges(94312388,244006886))

require(Homo.sapiens)
gns = geneRanges(Homo.sapiens, column="SYMBOL")

output = CNVanno(gns, cnv, "SYMBOL")

cnvMatanno<-output$cnvMatanno
cnvGeneMat<-output$cnvGeneMat

write.table(cnvMatanno,paste0(spl_name,"_cnvMatanno.txt"),sep="\t",row.names=F)
write.table(cnvGeneMat,paste0(spl_name,"_cnvGeneMat.txt"),quote=F,sep="\t",row.names=F)

print (pasteo("Copy number annotation of ",spl_name," Finshed!")
