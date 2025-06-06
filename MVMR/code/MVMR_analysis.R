######################################################
## MVMR analysis
## exposure: sleep and BMI
## step-1: match sleep and BMI
## step-2: match with outcome
## step-3: do harmonization
## step3-1: harmonize between exposures
## step3-2: harmonize between exposure and outcome
## step-4: MVMR analysis 
## note: for sleep and outcome
## the association is "without" adjusting for BMI
######################################################

library(TwoSampleMR)

library(data.table)

library(dplyr)


## load BMI summary statistics
BMI <- fread("/bmi.giant-ukbb.meta-analysis.females.23May2018.txt")

BMI$rsID <- sub(":.*", "", BMI$SNP)

head(BMI)


## exposure
exposure <- "OSA"
#exposure <- "Insomnia"
#exposure <- "ShortSleep"
#exposure <- "LongSleep"
#exposure <- "Sleepiness"


## load sleep GWAS
## without BMI adjustment
sleep <- fread(paste0("/MVMR/Sleep/",
                      exposure, "/", 
                      exposure, "_bmi_unadj_XXXXX_pvalue_10_5.csv"))

match_sleep_BMI <- left_join(sleep, BMI, by = "rsID")

## remove variants
## if the variant cannot be matched to BMI GWAS

match_sleep_BMI <- match_sleep_BMI[!is.na(match_sleep_BMI$SNP), ]

## outcome
outcome <- "AF"
#outcome <- "CAD"
#outcome <- "CKD"
#outcome <- "HF"
#outcome <- "HTN"
#outcome <- "T2DM"

## load outcome GWAS
## without BMI adjustment
Outcome <- fread(paste0("/MVMR/Outcome/",
                        outcome, "/", 
                        outcome, "_snp_association_XXXXX_bmi_unadj.csv"))

Outcome <- Outcome[which(Outcome$rsid %in% match_sleep_BMI$rsID), ]

match_sleep_BMI_outcome <- match_sleep_BMI[which(match_sleep_BMI$rsID %in% Outcome$rsid), ]

head(match_sleep_BMI_outcome)

######################################################
## important step
## first harmonize between exposures
## using harmonise_data() function
######################################################

## see exposure of interest as "exposure"
DH_exp <- match_sleep_BMI_outcome[, c("rsID", "OR", "logOR_se", "PValue", "A1", "A2", "EAF")]

## change to log OR scale
DH_exp$OR <- log(DH_exp$OR)

colnames(DH_exp) <- c("SNP", "beta.exposure", "se.exposure", "pvalue.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")

DH_exp$exposure <- exposure

## randomly set the id name
DH_exp$id.exposure <- 1

## see second exposures as "outcome"
DH_out <- match_sleep_BMI_outcome[,c("rsID", "BETA", "SE", "P", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele")]

colnames(DH_out) <- c("SNP", "beta.outcome", "se.outcome", "pvalue.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")

DH_out$outcome <- outcome

## randomly set the id name
DH_out$id.outcome <- 2

## do harmonization between exposures
exposure_DH <- harmonise_data(DH_exp, DH_out, action = 2)

## remove variants did not pass filtering
exposure_DH <- exposure_DH[which(exposure_DH$mr_keep == "TRUE"), ]

######################################################
## data preprocess 
## before doing harmonization
######################################################

write_gwas_file <- function(filename, snp, beta, se, eaf, effect_allele, other_allele, pvalue, trait, id, n){
  data <- data.frame(
    snp = snp,
    beta = beta,
    se = se,
    eaf = eaf,
    effect_allele = effect_allele,
    other_allele = other_allele,
    pvalue = pvalue,
    trait = trait,
    id = id,
    n = n)
  write.table(data, file = filename, quote = FALSE, row.names = FALSE, sep = " ")
}


## first exposure: sleep exposure
write_gwas_file(filename = paste0("exposure_", exposure, ".txt"), 
                snp = exposure_DH$SNP, 
                beta = exposure_DH$beta.exposure,
                se = exposure_DH$se.exposure,
                eaf = exposure_DH$eaf.exposure,
                effect_allele = exposure_DH$effect_allele.exposure,
                other_allele = exposure_DH$other_allele.exposure,
                pvalue = exposure_DH$pvalue.exposure,
                trait = exposure, id = paste0(exposure, "_gwas"), n = 100000)

## second exposure: BMI
write_gwas_file(filename = "exposure_BMI.txt", 
                snp = exposure_DH$SNP, 
                beta = exposure_DH$beta.outcome,
                se = exposure_DH$se.outcome,
                eaf = exposure_DH$eaf.outcome,
                effect_allele = exposure_DH$effect_allele.outcome,
                other_allele = exposure_DH$other_allele.outcome,
                pvalue = exposure_DH$pvalue.outcome,
                trait = "BMI", id = "BMI_gwas", n = 100000)

## outcome: CVD-outcome
write_gwas_file(filename = paste0("outcome_", outcome, ".txt"), 
                snp = Outcome$rsid,
                beta = Outcome$beta,
                se = Outcome$se,
                eaf = Outcome$female_AF_alt,
                effect_allele = Outcome$alt,
                other_allele = Outcome$ref,
                pvalue = Outcome$pvalue,
                trait = outcome, id = paste0(outcome, "_gwas"), n = Outcome$number_case)

######################################################
## load in exposure and outcome data
## in a format from TwoSampleMR package
## use "read_exposure_data" and "read_outcome_data" function
######################################################
######################################################
## load exposures
## use "read_exposure_data()" function
######################################################
exposure_files <- c(paste0("exposure_", exposure, ".txt"), "exposure_BMI.txt")
exposures <- lapply(exposure_files, function(XX){
  read_exposure_data(
    filename = XX,
    sep = " ",
    phenotype_col = "trait",
    snp_col = "snp",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pvalue",
    samplesize_col = "n",
    id_col = "id")}
)

exp_combined <- do.call(rbind, exposures)

######################################################
## load outcome
## use "read_outcome_data()" function
######################################################
outcome_data <- read_outcome_data(
  filename = paste0("outcome_", outcome, ".txt"),
  sep = " ",
  phenotype_col = "trait",
  snp_col = "snp",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pvalue",
  samplesize_col = "n",
  id_col = "id")


######################################################
## step-3
## harmonize exposure and outcome data
## use "mv_harmonise_data()" function
######################################################

MVMR_HM_data <- mv_harmonise_data(exposure_dat = exp_combined, outcome_dat = outcome_data)

######################################################
## step 4-1
## prepare the final data for MVMR analysis
## use "mr_mvinput()" function from "MendelianRandomization" R package
library(MendelianRandomization)
######################################################

Input_data <- mr_mvinput(bx = MVMR_HM_data$exposure_beta, 
                         bxse = MVMR_HM_data$exposure_se,
                         by = MVMR_HM_data$outcome_beta, 
                         byse = MVMR_HM_data$outcome_se)

######################################################
## step 4-2
## run MVMR
## use the method from "MendelianRandomization" R package
## IVW, robust IVW, median
######################################################

## Multivariable inverse-variance weighted method
MVIVW <- mr_mvivw(Input_data)

## robust Multivariable inverse-variance weighted method
## using robust regression
MVIVW_robust <- mr_mvivw(Input_data, robust = TRUE)

## Multivariable median-based method
MV_median <- mr_mvmedian(Input_data, iterations = 10000)

## save results
## MVMR-IVW
MR_MVIVW_res <- data.frame(MVIVW$class[1],
                           MVIVW$Estimate,
                           MVIVW$StdError,
                           MVIVW$CILower,
                           MVIVW$CIUpper,
                           MVIVW$Pvalue,
                           MVIVW$SNPs)

colnames(MR_MVIVW_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue", "number_IVs")


## save results
## MVMR-robust-IVW
MR_MVIVW_robust_res <- data.frame(paste0(MVIVW_robust$class[1], "_robust"),
                                  MVIVW_robust$Estimate,
                                  MVIVW_robust$StdError,
                                  MVIVW_robust$CILower,
                                  MVIVW_robust$CIUpper,
                                  MVIVW_robust$Pvalue,
                                  MVIVW_robust$SNPs)

colnames(MR_MVIVW_robust_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue", "number_IVs")


## save results
## MVMR-median
MV_median_res <- data.frame(MV_median$class[1],
                            MV_median$Estimate,
                            MV_median$StdError,
                            MV_median$CILower,
                            MV_median$CIUpper,
                            MV_median$Pvalue,
                            MV_median$SNPs)

colnames(MV_median_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue", "number_IVs")

## combine all results
result <- rbind(MR_MVIVW_res, MR_MVIVW_robust_res, MV_median_res)

######################################################
## step 5
## for checking IV strength
## use R package "MVMR"
######################################################

library(remotes)

#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

library(MVMR)

check_data <- matrix(0, length(MVMR_HM_data$outcome_beta), 7)

check_data[, c(1,2)] <- MVMR_HM_data$exposure_beta
check_data[, c(3,4)] <- MVMR_HM_data$exposure_se
check_data[, 5] <- MVMR_HM_data$outcome_beta
check_data[, 6] <- MVMR_HM_data$outcome_se
check_data[, 7] <- row.names(MVMR_HM_data$exposure_beta)

F.data <- format_mvmr(
  BXGs = check_data[, c(1, 2)],
  BYG = check_data[, 5],
  seBXGs = check_data[, c(3, 4)],
  seBYG = check_data[, 6],
  RSID = check_data[, 7])

## IV strength
sres <- strength_mvmr(r_input = F.data, gencov = 0)

F.data$exposure_X1 <- "BMI"
F.data$exposure_X2 <- exposure
F.data$Y <- outcome
F.data$F_statistic_X1 <- sres$exposure1
F.data$F_statistic_X2 <- sres$exposure2


colnames(F.data) <- c("SNP", "beta_outcome", "se_outcome",
                      "beta_exp1", "beta_exp2", 
                      "se_exp1", "se_exp2",
                      "exp1", "exp2", "outcome", 
                      "cond_F_stat_exp1", "cond_F_stat_exp2")


