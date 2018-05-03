#TCGA Sample Classification 2018-05-02
#TCGA assign the tumor samples with 0-9, and normal samples with 10-19
setwd("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/")
tcgamate<-as.data.frame(read.csv("./metadata/TCGA_BRCA_Clinical_Info.csv",h=T))#TCGA clinial info: RACE SAMPLE name SAMPLE ID
id<-read.table("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/TCGA_SNPSample_Annotation.txt",h=T)#sample name + TCGA ID

#PCAsampleinfo<-read.csv("clinical_PANCAN_patient_with_followup.tsv",h=T,sep="\t")
#table(PCAsampleinfo[,3])


SaplAnno<-cbind(gsub(".birdseed.data.txt","",id[,2]),substr(id[,7],1,12),substr(id[,7],14,16),c("Tumor","Normal")[as.factor(as.numeric(substr(id[,7],14,15))>=10)])
write.csv(SaplAnno,"TCGA_Array6_Sample_information.csv")

tumorsample<-SaplAnno[which(SaplAnno[,4]=="Tumor"),]
SumTumor<-as.matrix(table(tumorsample[,2]))
SumTumorSort<-SumTumor[order(SumTumor[,1],decreasing=T),]
tumorsample[grep("TCGA-A7-A26E",tumorsample[,2]),]


write.table(paste0(tumorsample[,1],".SNP6"),"tumorsamplename.txt",quote=F,row.names=F,col.names=F)
normalsample<-SaplAnno[which(SaplAnno[,4]=="Normal"),]
write.table(paste0(normalsample[,1],".SNP6"),"normalsamplename.txt",quote=F,row.names=F,col.names=F)

#01A  01B  06A  10A  10B  10C  10D  11A  11B 
#Sampletype: 01A  01B  06A  10A  10B  10C  10D  11A  11B 
#TypeNumber:1086   20    7  973   27   10    1  122   17 
#Code and tumor type
#01: Primary Solid Tumor
#06: Metastatic
#10: Blood Derived Normal
#11: Solid Tissue Normal



############################Select population with 1 T and 1 N fresh frozen, if no fresh frozen, use FFPE###############
##################################################Tumor################################################################
#duplicated samples
tmp1<-as.matrix(table(tumorsample[,2]))
dupspl<-rownames(tmp1)[which(tmp1[,1]>1)]
tumorsampleDup<-tumorsample[which(tumorsample[,2]%in%dupspl),]
tumorsampleDupSort<-tumorsampleDup[order(tumorsampleDup[,2]),]
#reasons of duplication:1. test both fresh frozen sample and FFPE sample; 2. test both primary tumor and metastasis tumor; 3. one tumor sample (TCGA-B6-A1KC) test twice 
#So, just remove the duplicated FFPE samples and metastasis samples
#[1,] "CUSKS_p_TCGAb47_SNP_1N_GenomeWideSNP_6_H08_628250"          "TCGA-A7-A0DC" "01A" "Tumor"
#[2,] "CHILE_p_TCGASNP_230_231_FFPE_N_GenomeWideSNP_6_A10_1147968" "TCGA-A7-A0DC" "01B" "Tumor" remove
#[3,] "CHILE_p_TCGASNP_230_231_FFPE_N_GenomeWideSNP_6_A07_1147964" "TCGA-A7-A13G" "01B" "Tumor" remove
#[4,] "RANDS_p_TCGA_b109_SNP_1N_GenomeWideSNP_6_A01_771906"        "TCGA-A7-A13G" "01A" "Tumor" 
#[5,] "BEDEW_p_TCGA_FFPE_7_13_N_GenomeWideSNP_6_C11_1347436"       "TCGA-A7-A26E" "01A" "Tumor" remove
#[6,] "BEDEW_p_TCGA_FFPE_7_13_N_GenomeWideSNP_6_C12_1347490"       "TCGA-A7-A26E" "01B" "Tumor" remove
#[7,] "EXINE_p_TCGA_b142_SNP_N_GenomeWideSNP_6_C12_802410"         "TCGA-A7-A26E" "01A" "Tumor"
#[8,] "CHILE_p_TCGASNP_230_231_FFPE_N_GenomeWideSNP_6_A09_1147944" "TCGA-A7-A26F" "01B" "Tumor" remove
#[9,] "EXINE_p_TCGA_b142_SNP_N_GenomeWideSNP_6_A12_802324"         "TCGA-A7-A26F" "01A" "Tumor"
#[10,] "BEDEW_p_TCGA_FFPE_7_13_N_GenomeWideSNP_6_C06_1347318"       "TCGA-A7-A26J" "01B" "Tumor" remove
#[11,] "EXINE_p_TCGA_b142_SNP_N_GenomeWideSNP_6_A09_802424"         "TCGA-A7-A26J" "01A" "Tumor"
#[12,] "BEDEW_p_TCGA_FFPE_7_13_N_GenomeWideSNP_6_C05_1347438"       "TCGA-A7-A26J" "01A" "Tumor" remove
#[13,] "CELLA_p_TCGA_177_179_Esph_SNP_N_GenomeWideSNP_6_E06_871070" "TCGA-AC-A2QH" "01A" "Tumor"
#[14,] "CHILE_p_TCGASNP_230_231_FFPE_N_GenomeWideSNP_6_G04_1147836" "TCGA-AC-A2QH" "01B" "Tumor" remove
#[15,] "CENTS_p_TCGASNP_212_216_217_N_GenomeWideSNP_6_C04_1039446"  "TCGA-AC-A3OD" "01A" "Tumor"
#[16,] "CHILE_p_TCGASNP_230_231_FFPE_N_GenomeWideSNP_6_E04_1147978" "TCGA-AC-A3OD" "01B" "Tumor" remove
#[17,] "FENCE_p_TCGAb_332_333_334_NSP_GenomeWideSNP_6_E09_1367264"  "TCGA-AC-A6IX" "01A" "Tumor" 
#[18,] "FENCE_p_TCGAb_332_333_334_NSP_GenomeWideSNP_6_E10_1367332"  "TCGA-AC-A6IX" "06A" "Tumor" Duplicate metastasis
#[19,] "RANDS_p_TCGA_b109_SNP_1N_GenomeWideSNP_6_A07_771954"        "TCGA-B6-A1KC" "01A" "Tumor"
#[20,] "USAGE_p_TCGA_132_133_136_SNP_N_GenomeWideSNP_6_F05_787566"  "TCGA-B6-A1KC" "01B" "Tumor" Duplicate
#[21,] "HITCH_p_TCGASNP_b93_N_GenomeWideSNP_6_G08_741502"           "TCGA-BH-A18V" "01A" "Tumor"
#[22,] "MINAE_p_TCGA_200_202_203_SNP_N_GenomeWideSNP_6_D05_955116"  "TCGA-BH-A18V" "06A" "Tumor" Duplicate metastasis
#[23,] "GAMMA_p_TCGA_b103_104_SNP_N_GenomeWideSNP_6_B08_755844"     "TCGA-BH-A1ES" "01A" "Tumor"
#[24,] "LEONE_p_TCGA_103_243_257_N_GenomeWideSNP_6_A01_1300650"     "TCGA-BH-A1ES" "06A" "Tumor" Duplicate metastasis
#[25,] "RANDS_p_TCGA_b109_SNP_1N_GenomeWideSNP_6_F05_771958"        "TCGA-BH-A1FE" "01A" "Tumor"
#[26,] "MINAE_p_TCGA_200_202_203_SNP_N_GenomeWideSNP_6_D04_955122"  "TCGA-BH-A1FE" "06A" "Tumor" Duplicate metastasis
#[27,] "HITCH_p_TCGASNP_b93_N_GenomeWideSNP_6_D08_741406"           "TCGA-E2-A15A" "01A" "Tumor"
#[28,] "HITCH_p_TCGASNP_b93_N_GenomeWideSNP_6_A07_741428"           "TCGA-E2-A15A" "06A" "Tumor" Duplicate metastasis
#[29,] "HITCH_p_TCGASNP_b93_N_GenomeWideSNP_6_B10_741426"           "TCGA-E2-A15E" "06A" "Tumor" Duplicate metastasis
#[30,] "HITCH_p_TCGASNP_b93_N_GenomeWideSNP_6_E02_741534"           "TCGA-E2-A15E" "01A" "Tumor"
#[31,] "JOUAL_p_TCGA_b96_SNP_N_GenomeWideSNP_6_A11_748018"          "TCGA-E2-A15K" "06A" "Tumor" Duplicate metastasis
#[32,] "JOUAL_p_TCGA_b96_SNP_N_GenomeWideSNP_6_C10_748046"          "TCGA-E2-A15K" "01A" "Tumor"

#extract the duplicated samples
#duplicate FFPE samples
ffpedup<-tumorsampleDupSort[grep("_FFPE_",tumorsampleDupSort[,1]),1]
#duplicate metastasis samples
metastasisdup<-tumorsampleDupSort[which(tumorsampleDupSort[,3]=="06A"),1]
#duplicate fresh frezon samples
FFdup<-"USAGE_p_TCGA_132_133_136_SNP_N_GenomeWideSNP_6_F05_787566"

allDup<-c(ffpedup,metastasisdup,FFdup)

tumorremoveDup<-tumorsample[-pmatch(allDup,tumorsample[,1]),]

##################################################Normal################################################################
#duplicated samples
tmp2<-as.matrix(table(normalsample[,2]))
normaldupspl<-rownames(tmp2)[which(tmp2[,1]>1)]
normalsampleDup<-normalsample[which(normalsample[,2]%in%normaldupspl),]
normalsampleDupSort<-normalsampleDup[order(normalsampleDup[,2]),]
#summarize normalsampleDupSort into a file duplicated_normal_samples.csv
normalDup<-read.csv("duplicated_normal_samples.csv",h=T)
#10A 11A 11B 
#58  51   1 
#reasons of duplication:1. test both blood sample (10A) and para normal tissue sample (11A); 2. sequenced both FFPE and FF para sample ; 3. one blood sample (TCGA-XX-A899) test twice 
#total 55 samples with duplicate, among them: 2 sample with Duplicate FFPE; 1 sample with Duplicated blood; 52 sample with duplicated para
#So, just remove the duplicated FFPE samples, blood and para samples
normalDupSpl<-normalDup[grep("Duplicate",normalDup[,6]),2]

normalremoveDup<-normalsample[-pmatch(normalDupSpl,normalsample[,1]),]

#keep the sample in both tumor and normal samples
overlap<-intersect(tumorremoveDup[,2],normalremoveDup[,2])


tumorremoveFin<-tumorremoveDup[which(tumorremoveDup[,2]%in%overlap),]
normalremoveFin<-normalremoveDup[which(normalremoveDup[,2]%in%overlap),]
colnames(normalremoveFin)<-c("filename","sampleID","tissueType","Type")
colnames(tumorremoveFin)<-c("filename","sampleID","tissueType","Type")
write.table(normalremoveFin,"normalsample_information.txt",row.names = F,quote = F)
write.table(tumorremoveFin,"tumorsample_information.txt",row.names = F,quote = F)


#read Europen sample name (base on PCA classification)
european<-unique(read.table("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/SNP/PopulationClustering/Europen_SampleID_in_TCGA_BRCA.txt",h=T)[,1])

EuropeanNormal<-normalremoveFin[which(normalremoveFin[,2]%in%european),]
EuropeanTumor<-tumorremoveFin[which(tumorremoveFin[,2]%in%european),]
europeanoverlap<-intersect(unique(EuropeanNormal[,2]),EuropeanTumor[,2])

EuropeanNormalFin<-EuropeanNormal[which(EuropeanNormal[,2]%in%europeanoverlap),]
EuropeanTumorFin<-EuropeanTumor[which(EuropeanTumor[,2]%in%europeanoverlap),]

write.table(EuropeanNormalFin,"EuropeanNormalSample_information.txt",row.names = F,quote = F,sep="\t")
write.table(EuropeanTumorFin,"EuropeanTumorSample_information.txt",row.names = F,quote = F,sep="\t")
