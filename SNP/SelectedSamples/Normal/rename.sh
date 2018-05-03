#extract sample for following analysis
awk '{print $1".SNP6\t"$1".SNP6"}' normalsample_information.txt > exracted_sample_for_following_analysis.txt
plink1.9 --file /home/tao/taoqing/projects/sgi/TCGAArray6SNP/normalsample/merge_tcga_brca_normal --keep  exracted_sample_for_following_analysis.txt --make-bed --out normal_sample_matched2018may03

#sample rename
awk '{print $1".SNP6\t"$1".SNP6\t"$2"\t"$2}' normalsample_information.txt > Normal_rename_recode.txt
plink1.9 --bfile normal_sample_matched2018may03 --update-ids Normal_rename_recode.txt --make-bed --out normal_sample_matched2018may03_rename



#extract sample for following analysis
awk '{print $1".SNP6\t"$1".SNP6"}' tumorsample_information.txt > exracted_sample_for_following_analysisi_tumor.txt
plink1.9 --file /home/tao/taoqing/projects/sgi/TCGAArray6SNP/merge_tcga_brca_2663_snp_plink --keep  exracted_sample_for_following_analysis_tumor.txt --make-bed --out tumor_sample_matched2018may03

#sample rename
awk '{print $1".SNP6\t"$1".SNP6\t"$2"\t"$2}' tumorsample_information.txt > Tumor_rename_recode.txt
plink1.9 --bfile tumor_sample_matched2018may03 --update-ids Tumor_rename_recode.txt --make-bed --out tumor_sample_matched2018may03_rename





