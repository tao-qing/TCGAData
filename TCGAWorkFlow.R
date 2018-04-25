#source('https://bioconductor.org/biocLite.R')
#biocLite('GenomicDataCommons')
library(GenomicDataCommons)

setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/")

#
manifest <- read.table("gdc_manifest_20180328_144107.txt",h=T)
file_uuids <- manifest$id

TCGAtranslateID = function(file_ids, legacy = TRUE) {
  info = files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

res = TCGAtranslateID(as.character(file_uuids))

IDanno<-cbind(manifest,res[pmatch(file_uuids,rownames(res)),])

write.table(IDanno,"TCGA_SampleName_Annotation.txt")

#totaly 7,837 sample
#SNP birdseed 2263
#CNV ismpolish 2263
#EXP rsem 2430 (Normalized/rsem RNASeq)
#EXP quantification 881 (RNA-seq)




###########Data April 13, 2018#####annotate_snp_sample_info###################
snpsample<-as.data.frame(read.table("metaData/snp_fileinfo.txt",h=F))

res = TCGAtranslateID(as.character(snpsample[,1]))

IDanno<-cbind(snpsample,res[pmatch(as.character(snpsample[,1]),rownames(res)),])

write.table(IDanno,"TCGA_SNPSample_Annotation.txt")

















