#summarize for the germline copy number variation
setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/CNV/Normal/")

merge_cnv<-read.table("merge_cnv_fin.txt",h=T)
#210917 cnvs

#filter cnv log2Ratio > 0.5 or < -0.5; Ratio > 1.4 or < 0.7
merge_cnv_filt<-merge_cnv[c(which(merge_cnv$log2Ratio>0.5),which(merge_cnv$log2Ratio<(-0.5))),]
#18575 cnvs

samples<-c(unique(merge_cnv$sample),"TCGA-BH-A0DL_10A","TCGA-D8-A27W_10A")
#1059 samples

genes<-unique(merge_cnv$Genes)
#19295 genes

library(foreach)
library(doParallel)
registerDoParallel(cores=100)

r1=foreach(ge=genes,.combine = "rbind")%dopar%{
  mat_loss<-merge_cnv[which(merge_cnv$cnvStatus=="Loss"),]
  mat_gain<-merge_cnv[which(merge_cnv$cnvStatus=="Gain"),]
  loss_freq<-length(unique(mat_loss$sample[grep(ge,mat_loss$Genes)]))/1059
  gain_freq<-length(unique(mat_gain$sample[grep(ge,mat_gain$Genes)]))/1059
  c(as.character(ge),loss_freq,gain_freq)
}

write.csv(r1,"gene_cnv_freq_filter0.5.csv")

cnvMat<-matrix(data=0,ncol=length(unique(merge_cnv$sample)),nrow=length(unique(merge_cnv$Genes)),dimnames=list(unique(merge_cnv$Genes),unique(merge_cnv$sample)))

i=1
for(ge in genes){
  index<-grep(ge,merge_cnv$Genes)
  spl<-merge_cnv$sample[index]
  ratio<-merge_cnv$log2Ratio[index]
  cnvMat[grep(ge,rownames(cnvMat)),sapply(spl,function(x)grep(x,colnames(cnvMat)))]<-ratio
  i=i+1
  print(i)
}

write.csv(cnvMat,"copynumbervariationMat_filter0.5.csv")

foreach(ge=genes,.combine = "rbind")%dopar%{
  index<-grep(ge,merge_cnv$Genes)
  spl<-merge_cnv$sample[index]
  ratio<-merge_cnv$log2Ratio[index]
  cnvMat[grep(ge,rownames(cnvMat)),sapply(spl,function(x)grep(x,colnames(cnvMat)))]<-ratio
  print(ge)
}


cnv_frep<-read.csv("gene_cnv_freq_filter0.5.csv",fill=NA)
cnv_frep<-cnv_frep[,-1]
cnv_frepSortLoss<-cnv_frep[order(cnv_frep[,2],decreasing=T),]

cosmicgene<-unique(read.csv("CosmicMutantExportCensus.tsv",sep="\t")[,1])


#DNA Damage Response Gene
ddrgene<-c("PARP1","POLB","POLG","POLD1","POLD1","POLE1","POLE","PCNA","REV3L","POLZ","MAD2L2","REV7","REV1L","REV1","POLH","POLI","RAD30B","POLQ","POLK","DINB1","POLL","POLM","POLN","POL4P")
#MMR genes
mmrgene<-c("MSH2","MSH3","MSH6","MLH1","PMS2","MSH4","MSH5","MLH3","PMS1","PMS2L3","SETD2","MYH11","EPCAM","TGFBR2", "MUTYH")
#BRCA pathway
brcagene<-c("FANCE","FANCC","FANCA","FANCL","FANCG","FANCF","FANCD2","ATM","ATR","CHK2","BRCA1","BRCA2","BARD1","BRCD1","BACH1","BLM","NBS1","MRE11","RAD50","SMARCA1","ATF1","CHK1","PLK1","CDKN1A","GADD45A","OTC1","RAD51")
#NER
NER<-c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8","DDB1","DDB2","GTF2H1","GTF2H2","GTF2H3","GTF2H4","GTF2H5","LIG1","RAD23A","RAD23B","XPA","XPC","CETN2","CUL4B","CUL4A","CDK7","MNAT1","UVSSA","MMS19","ERCC6","PGBD3","BIVM","ERCC5")

cancergene<-unique(c(ddrgene,mmrgene,brcagene,NER))

cancergenecnv<-cnv_frepSortLoss[which(cnv_frepSortLoss[,1]%in%cancergene),]


#################################Mutation Matrix######################################
cnvgeneMat<-round(read.csv("copynumbervariationMat_filter0.5.csv",h=T,r=1),digits=2)
library(ComplexHeatmap)
library(circlize)
#Sample=as.character(group)
#ha <- HeatmapAnnotation(as.data.frame(Sample), col = list(Sample=c("N" = "red","T" = "green")))
png("cnv_allgene_heatmap.png",height=1200,width=1500,res=96)
Heatmap(as.matrix(cnvgeneMat),col = colorRamp2(c(-0.5,0, 0.5), c( "green","floralwhite","red")),cluster_rows = TRUE, cluster_columns = TRUE,show_column_names=FALSE,show_row_names = FALSE,heatmap_legend_param = list(title = "log2(RPM)", ncol = 1))
dev.off()


#########################################Cancer Gene Heatmap###########################
rasgene<-c("BRAF","FOS","JUN","ILK","KRAS","EGFR","VEGFR","ITGB4","NF1","NF2","VHL","FH","SDHB")
cellcyclegene<-c("ABL1","ABL2","CDK2","CDK4","CCND1","CCNE1","CDKN2B","CDKN2A","CDKN1C","RB1","AURKA","MDM2","CDKN2A","ATM","ATR","BRCA1","CHEK1","CHEK2")
#onco gene
oncogene<-c("ABL1","ABL2","AKT1","AKT2","ATF1","BCL11A","BCL2","BCL3","BCL6","BCR","BRAF","CARD11","CBLB","CBLC","CCND1","CCND2","CCND3","CDX2","CTNNB1","DDB2","DDIT3","DDX6","DEK","EGFR","ELK4","ERBB2","ETV4","ETV6","EVI1","EWSR1","FEV","FGFR1","FGFR1OP","FGFR2","FUS","GOLGA5","GOPC","HMGA1","HMGA2","HRAS","IRF4","JUN","KIT","KRAS","LCK","LMO2","MAF","MAFB","MAML2","MDM2","MET","MITF","MLL","MPL","MYB","IDH1","MYC","MYCL1","MYCN","NCOA4","NFKB2","NRAS","NTRK1","NUP214","PAX8","PDGFB","PIK3CA","PIM1","PLAG1","PPARG","PTPN11","RAF1","REL","RET","ROS1","SMO","SS18","TCL1A","TET2","TFG","TLX1","TPR","USP6","TCF3")
tsgene<-c("APC","ARHGEF12","ATM","BCL11B","BLM","BMPR1A","BRCA1","BRCA2","CARS","CBFA2T3","CDH1","CDH11","CDK6","CDKN2C","CEBPA","CHEK2","CREB1","CREBBP","CYLD","DDX5","EXT1","EXT2","FBXW7","FH","FLT3","FOXP1","GPC3","IL2","JAK2","MAP2K4","MDM4","MEN1","MLH1","MSH2","NF1","NF2","NOTCH1","NPM1","NR4A3","NUP98","PALB2","PML","PTEN","RB1","RUNX1","SDHB","SDHD","SMARCA4","SMARCB1","SOCS1","STK11","SUFU","SUZ12","SYK","TNFAIP3","TP53","TSC1","TSC2","VHL","WRN","WT1")

brcagenes<-c("ACVRL1","AFF2","AGMO","AGTR2","AHNAK","AHNAK2","AKAP9","AKT1","AKT2","ALK","APC","ARID1A","ARID1B","ARID2","ARID5B","ASXL1","ASXL2","ATR","BAP1","BCAS3","BIRC6","BRAF","BRCA1","BRCA2","BRIP1","CACNA2D3","CASP8","CBFB","CCND3","CDH1","CDKN1B","CDKN2A","CHD1","CHEK2","CLK3","CLRN2","COL12A1","COL22A1","COL6A3","CTCF","CTNNA1","CTNNA3","DCAF4L2","DNAH11","DNAH2","DNAH5","DTWD2","EGFR","EP300","ERBB2","ERBB3","ERBB4","FAM20C","FANCA","FANCD2","FBXW7","FLT3","FOXO1","FOXO3","FOXP1","FRMD3","GATA3","GH1","GLDC","GPR124","GPR32","GPS2","HDAC9","HERC2","HIST1H2BC","HRAS","JAK1","KDM3A","KDM6A","KLRG1","KMT2C","KRAS","L1CAM","LAMA2","LAMB3","LARGE","LDLRAP1","LIFR","LIPI","MAGEA8","MAP2K4","MAP3K1","MAP3K10","MAP3K13","MBL2","MEN1","MLL2","MLLT4","MTAP","MUC16","MYH9","MYO1A","MYO3A","NCOA3","NCOR1","NCOR2","NDFIP1","NEK1","NF1","NF2","NOTCH1","NPNT","NR2F1","NR3C1","NRAS","NRG3","NT5E","OR6A2","PALLD","PBRM1","PDE4DIP","PIK3CA","PIK3R1","PPP2CB","PPP2R2A","PRKACG","PRKCE","PRKCQ","PRKCZ","PRKG1","PRPS2","PRR16","PTEN","PTPN22","PTPRD","PTPRM","RASGEF1B","RB1","ROS1","RPGR","RUNX1","RYR2","SBNO1","SETD1A","SETD2","SETDB1","SF3B1","SGCD","SHANK2","SIAH1","SIK1","SIK2","SMAD2","SMAD4","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SPACA1","STAB2","STK11","STMN2","SYNE1","TAF1","TAF4B","TBL1XR1","TBX3","TG","THADA","THSD7A","TP53","TTYH1","UBR5","USH2A","USP28","USP9X","UTRN","ZFP36L1")


geneCNVMap<-function(mat,label,genes){
  overlap<-intersect(genes,rownames(mat))
  cnvMatsub<-as.matrix(mat[pmatch(overlap,rownames(mat)),])
  cnvgeneMatSort<-cnvMatsub[,order(apply(cnvMatsub,2,function(x)sum(abs(x))),decreasing=T)]
  freq<-round((length(which(apply(cnvgeneMatSort,2,function(x)any(x!=0))))/1061)*100,digits=2)
  tiff(paste0("cnv_",label,"_heatmap.png"),height=500*length(genes)/100,width=1000,res=90)
  Heatmap(cnvgeneMatSort,col = colorRamp2(c(-0.5,0, 0.5), c( "green","floralwhite","red")),column_title=paste0(label," (",freq,"%)"),cluster_rows = TRUE, cluster_columns = FALSE,show_column_names=FALSE,show_row_names = TRUE,heatmap_legend_param = list(title = "log2(RPM)"))
  dev.off()
}


mat=cnvgeneMat;label="DNA damage response";genes=cancergene
mat=cnvgeneMat;label="Onco driver gene";genes=oncogene
mat=cnvgeneMat;label="Tumor supressor";genes=tsgene
mat=cnvgeneMat;label="Cell cycle";genes=cellcyclegene
mat=cnvgeneMat;label="Ras";genes=rasgene
mat=cnvgeneMat;label="COSMIC";genes=cosmicgene
mat=cnvgeneMat;label="Gene Associated with BRCA";genes=brcagenes


geneCNVMap(cnvgeneMat,"Onco driver gene",oncogene)
geneCNVMap(cnvgeneMat,"Tumor supressor",tsgene)
geneCNVMap(cnvgeneMat,"Cell cycle",cellcyclegene)
geneCNVMap(cnvgeneMat,label="Ras",genes=rasgene)






