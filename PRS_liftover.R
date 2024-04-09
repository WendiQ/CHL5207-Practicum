# liftover PRS weights

library(data.table)

prs_weights <- fread("/Users/wendi qu/Documents/Data/PRS/PRS_weights.txt")
lift_fail <- fread("/Users/wendi qu/Documents/Data/PRS/lift_to_exclude.txt")

prs_post_exclude <- anti_join(prs_weights, lift_fail, by = c("chr", "position_hg19"))

# merge in the liftover SNPs
lifted <- fread("/Users/wendi qu/Documents/Data/PRS/PRS_lifted.txt")

prs_post_exclude <- prs_post_exclude %>%
  rename(chr_hg19 = chr) 

# append the lifted SNPs
prs_lifted <- cbind(lifted, prs_post_exclude)

prs_lifted<- prs_lifted %>%
  mutate(SNP = paste0("chr",CHR,":", POS, ":", A1, ":", A2))

genotype <- fread("/Users/wendi qu/Documents/Data/Merged-TOPMed-7arrays-10XG-PacBioHIFI-N5467-n7217886-postQC.bim") # 7217886

common <- merge(prs_lifted, genotype, by.x = c("CHR", "POS"), by.y = c("V1", "V4"))
# some are different in terms of effect alleles

# add column to both that include only chr:pos
prs_lifted <- prs_lifted %>%
  mutate(SNP = paste0("chr",CHR,":", POS))

genotype <- genotype %>%
 mutate(V2 = paste0("chr",V1,":", V4)) 
# 7215869
dups <- genotype[duplicated(genotype$V2),]

common_weights <- common %>%
  select(V2, effect_allele, effect_weight)

write.table(prs_lifted, "/Users/wendi qu/Documents/Data/PRS/PRS_weights_hg38.txt",
            sep="\t", row.names = F, quote = FALSE)
write.table(common_weights, "/Users/wendi qu/Documents/Data/PRS/PRS_weights_hg38_common.txt",
            sep="\t", row.names = F, quote = FALSE)


########## load in PRS results #############
prs_score <- fread("/Users/wendi qu/Documents/Data/PRS/PRS_gen_cf.profile")
pheco <- fread("/Users/wendi qu/Documents/Data/phenotype_raw_bmi.txt")
prs_pheno <- merge(pheco, prs_score)

hist(prs_pheno$SCORESUM, xlab = "PRS")

bmi_prs <- lm(bmi ~ SCORESUM, data = prs_pheno)
summary(bmi_prs)
plot(prs_pheno$SCORESUM, prs_pheno$bmi)
abline(bmi_prs)

bmi_prs_cov <- lm(bmi ~ age + agesq + sex + PIPS, data = prs_pheno)
summary(bmi_prs_cov) # r2 = 0.1423

bmi_prs_cov1 <- lm(bmi ~ SCORESUM + age + agesq + sex + PIPS, data = prs_pheno)
summary(bmi_prs_cov1) # r2 = 0.1512
plot(prs_pheno$SCORESUM, prs_pheno$bmi, xlab = "PRS", ylab = "BMI in CF")
abline(bmi_prs_cov1)

bmi_prs_cov2 <- lm(bmi ~ age + agesq + sex, data = prs_pheno)
summary(bmi_prs_cov2) # r2 = 0.1032

bmi_prs_cov3 <- lm(bmi ~ SCORESUM + age + agesq + sex, data = prs_pheno)
summary(bmi_prs_cov3) # r2 = 0.1133

bmi_prs2 <- lm(bmi ~ SCORE, data = prs_pheno2)
summary(bmi_prs2)
plot(prs_pheno2$SCORE, prs_pheno2$bmi)
abline(bmi_prs2)

hist(prs_pheno$SCORESUM)

