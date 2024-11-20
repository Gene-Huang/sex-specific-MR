##############################################################
## TwoSampleMR data analysis
## the results of male-specific causal effect estimate
## exposure GWAS: MVP sleep phenotype 
## outcome GWAS: AoU CVD-related phenotype
## note: for each exposure-outcome combination, we have
## (i) female-specific and male-specific GWAS, and
## (ii) GWAS with bmi-adjusted and bmi-unadjusted 
##############################################################

## load R packages
library(data.table)

library(dplyr)

## load MR function
source("/TwoSampleMR_methods_function.R")

## load SEB function
source("/SEB_function.R")

## load AW function
source("/Freq_AW_function.R")

## get variable (arguments) name
args <- commandArgs(trailingOnly=TRUE)

## OSA, Insomnia, Sleepiness, ShortSleep, LongSleep
exposure.pheno <- args[1]

## bmi.type = "bmi_unadj", "bmi_adj"
bmi.type <- args[2]

## AF, CAD, CKD, HF, HTN, T2DM
Outcome <- args[3]

## exposure GWAS = "male", "female"
sex <- "male"

## exposure GWAS pvalue = "5_10_8", "10_7", "10_5"
pvalue <- "10_5"

## get the sample size of (exposure) GWAS
## the input of cML method
Input_cML <- read.csv("/exposure_GWAS_sample_size.csv", header = T)

#####################################################
## start MR analysis
#####################################################
#####################################################
## step1: prepare input data
#####################################################

## read in harmonized exposure/outcome summary statistics
## select the final variants for MR analysis
harmonized_data <- fread("/Harmonized_exposure_outcome.csv")

## read in female exposure GWAS summary statistics
Exp.GWAS.F <- fread("/female_exposure_GWAS.csv")

Exp.GWAS.F <- Exp.GWAS.F[which(Exp.GWAS.F$SNP %in% harmonized_data$SNP), ]

## read in male exposure GWAS summary statistics
Exp.GWAS.M <- fread("/male_exposure_GWAS.csv")

Exp.GWAS.M <- Exp.GWAS.M[which(Exp.GWAS.M$SNP %in% harmonized_data$SNP), ]

## read in male outcome GWAS summary statistics
outcome.GWAS <- fread("male_exposure_GWAS.csv")

outcome.GWAS <- outcome.GWAS[which(outcome.GWAS$SNP %in% harmonized_data$SNP), ]

## unify the column name
colnames(Exp.GWAS.F) <- c("beta", "se", "pval")

## unify the column name
colnames(Exp.GWAS.M) <- c("beta", "se", "pval")

#####################################################
## step2: compute shrinkage estimate
## male-specific shrinkage estimate
## only apply for variant-exposure effect estimate
#####################################################

## meta
meta.M <- SEB(target_effect = Exp.GWAS.M$beta,
              other_effect = Exp.GWAS.F$beta,
              sd_target = Exp.GWAS.M$se,
              sd_other = Exp.GWAS.F$se,
              method = "meta")

## APM
APM.M <- SEB(target_effect = Exp.GWAS.M$beta,
             other_effect = Exp.GWAS.F$beta,
             sd_target = Exp.GWAS.M$se,
             sd_other = Exp.GWAS.F$se,
             method = "APM")

## freq AW
AW.M <- AW.Est.fn(target_effect = Exp.GWAS.M$beta,
                  other_effect = Exp.GWAS.F$beta,
                  sd_target = Exp.GWAS.M$se,
                  sd_other = Exp.GWAS.F$se)


#####################################################
## step3: apply twoSampleMR methods
## compute male-specific causal effect estimate
#####################################################
## variant-exposure: raw
res.M <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                             Outcome_se = outcome.GWAS$se,
                             Exposure_effect = Exp.GWAS.M$beta,
                             Exposure_se = Exp.GWAS.M$se,
                             N = Input_cML)


## variant-exposure: meta
res.meta.M <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                                  Outcome_se = outcome.GWAS$se,
                                  Exposure_effect = meta.M$Est,
                                  Exposure_se = meta.M$sd,
                                  N = Input_cML)

## variant-exposure: APM
res.APM.M <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                                 Outcome_se = outcome.GWAS$se,
                                 Exposure_effect = APM.M$Est,
                                 Exposure_se = APM.M$sd,
                                 N = Input_cML)


## variant-exposure: AW
res.AW.M <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                                Outcome_se = outcome.GWAS$se,
                                Exposure_effect = AW.M$AW.Est,
                                Exposure_se = AW.M$sd.AW.Est,
                                N = Input_cML)


#####################################################
## finalize results
#####################################################
Final.res <- rbind(res.M, res.meta.M, res.APM.M, res.AW.M)

## add one column: the input of exposure GWAS
## we consider 9 two sample MR methods
Final.res$exp_GWAS <- rep(c("raw", "meta", "APM", "AW"), each = 9)

## add one column: female or male results
Final.res$sex <- sex

## add one column: exposure variable
Final.res$exposure <- exposure.pheno

## add one column: outcome variable
Final.res$outcome <- Outcome

## add one column: exposure GWAS bmi-adj/unadj
Final.res$type <- bmi.type

#####################################################
## save results
#####################################################
## save the information
save.folder <- "/output/"

save.file <- paste0(save.folder,
                    exposure.pheno, "/",
                    bmi.type, "/",
                    Outcome, "/",
                    exposure.pheno, "_", Outcome, "_", sex, "_pvalue_", pvalue, ".csv")

write.csv(Final.res, save.file, row.names = F)

