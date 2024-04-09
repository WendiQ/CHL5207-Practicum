#!/bin/bash 


cd /hpf/largeprojects/struglis/wendi
module load plink/2.0
module load plink/1.90b6.21


####################################### Load in genotype files ########################################
plink2 --bfile /hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Combined-7Arrays-10XG-PacBioHiFi-PLINK/Merged-TOPMed-7arrays-10XG-PacBioHIFI-N5467-n7217886-postQC \
--keep subjects_keep.txt \
--not-chr x \
--extract SNP_keep.txt \
--make-bed --out cf_CAN # 2205


num_snps=$(wc -l < cf_CAN.bim)

# Print the result
echo "Number of SNPs in cf_CAN: $num_snps" # 7011436 47258


# Drop those with MAF < 0.05
plink --bfile cf_CAN --maf 0.05 --make-bed --out maf_CAN # 4472770 46230


####################################### Pruning high LD regions ########################################

plink --bfile maf_CAN \
--exclude High_LD_Regions_Build38.txt --range \
--make-bed --out LD_CAN # 4348502 45090

plink --bfile LD_CAN \
--indep-pairwise 1500 100 0.2 \
--out pruned_CAN_
plink --bfile LD_CAN \
--exclude pruned_CAN_.prune.out \
--make-bed --out pruned_CAN # 136795 32956

num_snps=$(wc -l < pruned_CAN.bim)

# Print the result
echo "Number of SNPs in pruned_CAN: $num_snps" # 32956


####################################### PCA ########################################

module load R/4.3.0
R


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GWASTools")
BiocManager::install("GENESIS")
BiocManager::install("gdsfmt", force=T)
BiocManager::install("SNPRelate", force=T)

library(GWASTools)
library(GENESIS)
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
setwd("/hpf/largeprojects/struglis/wendi")

# Convert PLINK binary files to GDS format
snpgdsBED2GDS(bed.fn = "pruned_CAN.bed",
bim.fn = "pruned_CAN.bim",
fam.fn = "pruned_CAN.fam",
out.gdsfn = "pruned_CAN.gds")

# Create a GenotypeData object from the GDS file
geno <- GdsGenotypeReader("pruned_CAN.gds")
genoData <- GenotypeData(geno)

# Perform principal component analysis (PCA)
mypcair <- pcair(genoData,verbose = TRUE)

# Extract eigenvalues and eigenvectors
eigenvalues <- mypcair$values
eigenvectors <- mypcair$vectors

# save eigenvalues and eigenvectors
write.table(eigenvalues, "PCAIR_CAN_Eigenvalues", quote=F, row.names=T, col.names=F)
write.table(eigenvectors, "PCAIR_CAN_Eigenvectors", quote=F, row.names=T, col.names=F)
write.table(mypcair$varprop[1:32], "PCAIR_CAN_varprop", quote=F, row.names=T, col.names=F)



####################### Previous code used for merging genotyped datasets #########################
###### H610 #########
prefix="/hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-H610/1_H610-N1722-flipped-all-chrs.chr"
suffix=".vcf.gz"
start_number=1
end_number=22

# Loop through the range of numbers
for ((number = start_number; number <= end_number; number++)); do
    # Construct the file name using the current number
    file_name="${prefix}${number}${suffix}"

    # Run PLINK command for each file
    plink --vcf "$file_name" --make-bed --out H610_chr${number}
done

plink2 --vcf /hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-H610/1_H610-N1722-flipped-all-chrs.chrX.vcf.gz \
--make-bed --out H610_chrX

# merge chr bfiles together
plink --merge-list H610_merge.txt --make-bed --out H610_merged

# Number of subjects included
plink --bfile H610_merged --hardy # 1722

###### H660 ########
prefix="/hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-H660/1_H660-N297-flipped-all-chrs.chr"
suffix=".vcf.gz"
start_number=1
end_number=22

# Loop through the range of numbers
for ((number = start_number; number <= end_number; number++)); do
    # Construct the file name using the current number
    file_name="${prefix}${number}${suffix}"

    # Run PLINK command for each file
    plink --vcf "$file_name" --make-bed --out H660_chr${number}
done

plink2 --vcf /hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-H660/1_H660-N297-flipped-all-chrs.chrX.vcf.gz \
--make-bed --out H660_chrX

# merge chr bfiles together
plink --merge-list H660_merge.txt --make-bed --out H660_merged

# Number of subjects included
plink --bfile H660_merged --hardy # 297

###### Omni5 ########
prefix="/hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-Omni5/2_Omni5-N61-flipped-all-chrs.chr"
suffix=".vcf.gz"
start_number=1
end_number=22

# Loop through the range of numbers
for ((number = start_number; number <= end_number; number++)); do
    # Construct the file name using the current number
    file_name="${prefix}${number}${suffix}"

    # Run PLINK command for each file
    plink --vcf "$file_name" --make-bed --out Omni5_chr${number}
done

plink2 --vcf /hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-Omni5/2_Omni5-N61-flipped-all-chrs.chrX.vcf.gz \
--make-bed --out Omni5_chrX

# merge chr bfiles together
plink --merge-list Omni5_merge.txt --make-bed --out Omni5_merged

# Number of subjects included
plink --bfile Omni5_merged --hardy # 61

###### Omni25 ########
prefix="/hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-Omni25/2_Omni25-N799-flipped-all-chrs.chr"
suffix=".vcf.gz"
start_number=1
end_number=22

# Loop through the range of numbers
for ((number = start_number; number <= end_number; number++)); do
    # Construct the file name using the current number
    file_name="${prefix}${number}${suffix}"

    # Run PLINK command for each file
    plink --vcf "$file_name" --make-bed --out Omni25_chr${number}
done

plink2 --vcf /hpf/largeprojects/struglis/CANFRE_consortium_GWAS_2023/Pre-TOPMed-processed-VCFs/TOPMed-Omni25/2_Omni25-N799-flipped-all-chrs.chrX.vcf.gz \
--make-bed --out Omni25_chrX

# merge chr bfiles together
plink --merge-list Omni25_merge.txt --make-bed --out Omni25_merged

# Number of subjects included
plink --bfile Omni25_merged --hardy # 799

rm H610_chr*
rm H660_chr*
rm Omni5_chr*
rm Omni25_chr*
rm H610_select*
rm H660_select*
rm Omni5_select*
rm Omni25_select*


# select SNPs that are in common for all arrays to merge them later
plink2 --bfile H610_merged \
--extract SNP_keep.txt \
--keep subjects_keep.txt \
--make-bed --out H610_select # 1312


plink2 --bfile H660_merged \
--extract SNP_keep.txt \
--keep subjects_keep.txt \
--make-bed --out H660_select # 164

plink2 --bfile Omni5_merged \
--extract SNP_keep.txt \
--keep subjects_keep.txt \
--make-bed --out Omni5_select # 38

plink2 --bfile Omni25_merged \
--extract SNP_keep.txt \
--keep subjects_keep.txt \
--make-bed --out Omni25_select # 457

plink --bfile H610_select \
--merge-list arrays.txt \
--make-bed --out merged_CAN1

# Exclude inconsistent SNPs and those with multiple positions

## Create a list of SNPs with multiple positions
grep 'Warning: Multiple positions seen for variant' merged_CAN1.log | awk '{print $NF}' | grep -o 'rs[0-9]*' > multiple_snps.txt

## combine inconsistent SNPs and multiple positions
cat multiple_snps.txt merged_CAN1-merge.missnp > snps_to_exclude.txt

plink --bfile H610_select --exclude snps_to_exclude.txt --make-bed --out H610_clean

plink --bfile H660_select --exclude snps_to_exclude.txt --make-bed --out H660_clean

plink --bfile Omni5_select --exclude snps_to_exclude.txt --make-bed --out Omni5_clean

plink --bfile Omni25_select --exclude snps_to_exclude.txt --make-bed --out Omni25_clean

plink --bfile H610_clean \
--merge-list clean_arrays.txt \
--make-bed --out merged_CAN





