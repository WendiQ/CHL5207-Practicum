library(data.table)
library(dplyr)

saige_result <- data.frame()

# combine output from saige
for (i in 1:22) {

  path <- paste0('/Users/wendi qu/Documents/Data/SAIGE chr/saige_step2_', i, '.txt')
  chr <- fread(path, header=T)
  saige_result <- rbind(saige_result, chr)
  
}

saige_result$chr <- as.numeric(gsub("\\D", "", saige_result$CHR))

saige_result <- fread("/Users/wendi qu/Documents/Data/SAIGE chr/saige_result.txt")

# add snp ID for locusfocus : e.g., 1_205720483_G_A_b37
saige_result <- saige_result %>%
  mutate(SNP = paste0(chr, "_", POS, "_", Allele1,"_",  Allele2, "_b38"))

saige_red <- saige_result %>%
  select(-c(MissingRate, Tstat, var, N, CHR, AC_Allele2, MarkerID))

write.table(saige_result, sep = "\t", file = "/Users/wendi qu/Documents/Data/SAIGE chr/saige_result.txt", 
            row.names = F, quote = FALSE)

write.table(saige_red, sep = "\t", file = "/Users/wendi qu/Documents/Data/SAIGE chr/saige_result_reduced.txt", 
            row.names = F, quote = FALSE)

# Manhattan plot 
manhattan(saige_result, chr="chr", bp="POS", snp="MarkerID", p="p.value" )


#################################### Sensitivity #######################################3

####### removing those with missing CFTR score #############
saige_result_sens <- data.frame()

# combine output from saige
for (i in 1:22) {
  
  path <- paste0('/Users/wendi qu/Documents/Data/Sensitivity/saige_step2_', i, '_nomiss.txt')
  chr <- read.table(path, header=T)
  saige_result_sens <- rbind(saige_result_sens, chr)
  
}

chr <- read.table('/Users/wendi qu/Documents/Data/Sensitivity/saige_step2_X_nomiss.txt', header=T)
saige_result_sens <- rbind(saige_result_sens, chr) # 5902465

saige_result_sens$chr <- as.numeric(gsub("\\D", "", saige_result_sens$CHR))
saige_result_sens$chr <- ifelse(is.na(saige_result_sens$chr), 23, saige_result_sens$chr) 

write.table(saige_result_sens, sep = "\t", file = "/Users/wendi qu/Documents/Data/Sensitivity/saige_result_sens.txt", 
            row.names = F, quote = FALSE)


############# Removing CFTR score #############
saige_result_noscore <- data.frame()

# combine output from saige
for (i in 1:22) {
  
  path <- paste0('/Users/wendi qu/Documents/Data/Sens_noscore/saige_step2_', i, '_noscore.txt')
  chr <- fread(path, header=T)
  saige_result_noscore <- rbind(saige_result_noscore, chr)
  
}


saige_result_noscore$chr <- as.numeric(gsub("\\D", "", saige_result_noscore$CHR))

write.table(saige_result_noscore, sep = "\t", file = "/Users/wendi qu/Documents/Data/Sens_noscore/saige_result_noscore.txt", 
            row.names = F, quote = FALSE)

############### including PI/PS status instead of CFTR score ##############
saige_result_pips <- data.frame()

# combine output from saige
for (i in 1:22) {
  
  path <- paste0('/Users/wendi qu/Documents/Data/Sens_pips/saige_step2_', i, '_pips.txt')
  chr <- read.table(path, header=T)
  saige_result_pips <- rbind(saige_result_pips, chr)
  
}


saige_result_pips$chr <- as.numeric(gsub("\\D", "", saige_result_pips$CHR))

# add snp ID for locusfocus : e.g., 1_205720483_G_A_b37
saige_pips_red <- saige_result_pips %>%
  mutate(SNP = paste0(chr, "_", POS, "_", Allele1,"_",  Allele2, "_b38")) %>%
  select(-c(MissingRate, Tstat, var, N, CHR, AC_Allele2, MarkerID))

write.table(saige_result_pips, sep = "\t", file = "/Users/wendi qu/Documents/Data/Sens_pips/saige_result_pips.txt", 
            row.names = F, quote = FALSE)
write.table(saige_pips_red, sep = "\t", file = "/Users/wendi qu/Documents/Data/Sens_pips/saige_result_pips_reduced.txt", 
            row.names = F, quote = FALSE)

saige_result_pips <- fread("/Users/wendi qu/Documents/Data/Sens_pips/saige_result_pips.txt")
############### including both PI/PS status and CFTR score ##############
saige_result_both_scores <- data.frame()

# combine output from saige
for (i in 1:22) {
  
  path <- paste0('/Users/wendi qu/Documents/Data/Sens_both_scores/saige_step2_', i, '_both.txt')
  chr <- fread(path, header=T)
  saige_result_both_scores <- rbind(saige_result_both_scores, chr)
  
}


saige_result_both_scores$chr <- as.numeric(gsub("\\D", "", saige_result_both_scores$CHR))

write.table(saige_result_both_scores, sep = "\t", file = "/Users/wendi qu/Documents/Data/Sens_both_scores/saige_result_both_scores.txt", 
            row.names = F, quote = FALSE)


############### CFTR score with NA values ##############
saige_result_nomiss <- data.frame()

# combine output from saige
for (i in 1:22) {
  
  path <- paste0('/Users/wendi qu/Documents/Data/Sens_nomiss/saige_step2_', i, '_nomiss.txt')
  chr <- read.table(path, header=T)
  saige_result_nomiss <- rbind(saige_result_nomiss, chr)
  
}

chr <- read.table('/Users/wendi qu/Documents/Data/Sens_nomiss/saige_step2_X_nomiss.txt', header=T)
saige_result_nomiss <- rbind(saige_result_nomiss, chr) # 5902465

saige_result_nomiss$chr <- as.numeric(gsub("\\D", "", saige_result_nomiss$CHR))
saige_result_nomiss$chr <- ifelse(is.na(saige_result_nomiss$chr), 23, saige_result_nomiss$chr)

write.table(saige_result_nomiss, sep = "\t", file = "/Users/wendi qu/Documents/Data/Sens_nomiss/saige_result_nomiss.txt", 
            row.names = F, quote = FALSE)



############################ Compare with general population ############################
saige_cf <- fread("/Users/wendi qu/Documents/Data/Sens_pips/saige_result_pips.txt", header=T)
gwas_gen <- fread("/Users/wendi qu/Documents/bmi.giant-ukbb.meta-analysis.combined.23May2018")
lift_fail <- fread("/Users/wendi qu/Documents/Data/lift_to_exclude.txt")
gen_build37_1 <- fread("/Users/wendi qu/Documents/Data/gen_build37.txt", header = F)
gen_build37_2 <- fread("/Users/wendi qu/Documents/Data/gen_build37_2.txt", header = F)

# remove SNPs in lift_fail from gwas_gen
gen_post_exclude <- anti_join(gwas_gen, lift_fail, by = c("CHR", "POS"))

# merge in the liftover SNPs
lifted <- fread("/Users/wendi qu/Documents/Data/lifted.txt")

gen_post_exclude_nox <- gen_post_exclude %>%
  rename(CHR37 = CHR,
         POS37 = POS) %>%
  filter(!(is.na(CHR37)))

# append the lifted SNPs
gen_lifted <- cbind(lifted, gen_post_exclude_nox)
write.table(gen_lifted, "/Users/wendi qu/Documents/Data/gen_lifted.txt")

## find SNPs in common
saige_cf$chr <- as.character(saige_cf$chr)
common_snps <- merge(saige_cf, gen_lifted, by.x = c("chr", "POS"), by.y = c("CHR", "POS")) # 4893919
write.table(common_snps, "/Users/wendi qu/Documents/Data/Common_snps_gen_cf.txt")

# check if snps significant in general population is also significant in CF
both_sig <- common_snps %>%
  filter(p.value < 0.05 & P < 5*10^(-8))
write.csv(both_sig, "/Users/wendi qu/Documents/Data/gen_cf_bothsig.csv")


cf_sig <- common_snps %>%
  filter(p.value < 10^(-5) & P >= 10^(-5))

gen_sig <- common_snps %>%
  filter(p.value >= 10^(-5) & P < 10^(-5))

sig_cf <- saige_cf[saige_cf$p.value < 10^(-5),]
sig_gen <- gwas_gen[gwas_gen$P < 10^(-5),]

gen_miss_sig <- sig_gen %>%
  filter()

# QQ plot
plot(-log10(common_snps$p.value), -log10(common_snps$P),
     xlab = "P-value in CF population", ylab = "P-value in general population")

# Draw a blue vertical line at x = 1.3
abline(v = 1.3, col = "blue")

# Draw a red vertical line at x = 7.3
abline(v = 7.3, col = "red")

# Draw a red horizontal line at y = 7.3
abline(h = 7.3, col = "red")

plot(-log10(common_snps$p.value), -log10(common_snps$P), ylim = c(0,100),
     xlab = "P-value in CF population", ylab = "P-value in general population")

# Draw a blue vertical line at x = 1.3
abline(v = 1.3, col = "blue")

# Draw a red vertical line at x = 7.3
abline(v = 7.3, col = "red")

# Draw a red horizontal line at y = 7.3
abline(h = 7.3, col = "red")


# x=1.3 (blue) x=7.3 (red) and a horizontal line at y=7.3 (red)
# check if top SNPs in CF are also sig in gen pop



