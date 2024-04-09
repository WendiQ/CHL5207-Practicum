#!/bin/bash

cd /hpf/largeprojects/struglis/wendi

module load R/4.2.2
module load htslib

Rscript step1_fitNULLGLMM.R     \
        --plinkFile=pruned_CAN  \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=pheno_covar.txt \
        --phenoCol=bmi \
        --covarColList=age,agesq,sex,PIPS,PC1,PC2,PC3,PC4,PC5,PC6 \
        --qCovarColList=sex  \
        --sampleIDColinphenoFile=IID \
        --invNormalize=FALSE     \
        --traitType=quantitative        \
        --outputPrefix=Sens_pips/saige_step1_pips \
        --nThreads=24	\
        --IsOverwriteVarianceRatioFile=TRUE

prefix="Dosage_data/Merged-PacBioHiFi-10XG-TOPMed7arrays-N5467-chr"
suffix=".vcf.gz"
suffix_index=".csi"
start_number=1
end_number=22

# Loop through the range of numbers
for ((number = start_number; number <= end_number; number++)); do
    # Construct the file name using the current number
    file_name="${prefix}${number}${suffix}"

    # Create index for each file
    tabix --csi -p vcf "$file_name" 

    # Run SAIGE step 2
    Rscript step2_SPAtests.R        \
        --vcfFile="$file_name"    \
        --vcfFileIndex="$file_name.csi" \
        --vcfField=DS   \
        --SAIGEOutputFile=Sens_pips/saige_step2_${number}_pips.txt    \
        --chrom=chr${number}       \
        --minMAF=0.01 \
        --GMMATmodelFile=Sens_pips/saige_step1_pips.rda \
        --varianceRatioFile=Sens_pips/saige_step1_pips.varianceRatio.txt \
        --is_output_moreDetails=TRUE    \
        --LOCO=TRUE

done



# tabix --csi -p vcf Dosage_data/Merged-PacBioHiFi-10XG-TOPMed7arrays-N5467-chrX-FINAL.vcf.gz 
# # zcat Merged-PacBioHiFi-10XG-TOPMed7arrays-N5467-chr21.vcf.gz  | less -S
# # Merged-PacBioHiFi-10XG-TOPMed7arrays-N5467-chr21.vcf.gz.csi


# Rscript step2_SPAtests.R        \
#         --vcfFile=Dosage_data/Merged-PacBioHiFi-10XG-TOPMed7arrays-N5467-chrX-FINAL.vcf.gz    \
#         --vcfFileIndex=Dosage_data/Merged-PacBioHiFi-10XG-TOPMed7arrays-N5467-chrX-FINAL.vcf.gz.csi \
#         --vcfField=DS   \
#         --SAIGEOutputFile=Sens_pips/saige_step2_X_pips.txt    \
#         --chrom=chrX       \
#         --minMAF=0.01 \
#         --GMMATmodelFile=Sens_pips/saige_step1_pips.rda \
#         --varianceRatioFile=Sens_pips/saige_step1_pips.varianceRatio.txt \
#         --is_output_moreDetails=TRUE    \
#         --LOCO=TRUE