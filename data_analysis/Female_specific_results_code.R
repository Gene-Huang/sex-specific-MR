##############################################################
## TwoSampleMR data analysis
## the results of female-specific causal effect estimate
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

## get variable name
args <- commandArgs(trailingOnly=TRUE)

## OSA, Insomnia, Sleepiness, ShortSleep, LongSleep
exposure.pheno <- args[1]

#bmi.type = "bmi_unadj", "bmi_adj"
bmi.type <- args[2]

## AF, CAD, CKD, HF, HTN, T2DM
Outcome <- args[3]

#exposure GWAS = "male", "female"
sex <- "female"

#exposure GWAS pvalue = "5_10_8", "10_7", "10_5"
pvalue <- "10_5"

## get the sample size of GWAS
## the input of cML method
Input_cML <- read.csv("/exposure_GWAS_sample_size.csv", header = T)

#####################################################
## start MR analysis
#####################################################
#####################################################
## prepare input data
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

## read in female outcome GWAS summary statistics
outcome.GWAS <- fread("female_exposure_GWAS.csv")

outcome.GWAS <- outcome.GWAS[which(outcome.GWAS$SNP %in% harmonized_data$SNP), ]

## unify the column name
colnames(Exp.GWAS.F) <- c("beta", "se", "pval")

## unify the column name
colnames(Exp.GWAS.M) <- c("beta", "se", "pval")

#####################################################
## compute female-specific shrinkage estimate
## only for variant-exposure effect estimate
#####################################################

## meta
meta.F <- SEB(target_effect = Exp.GWAS.F$beta,
              other_effect = Exp.GWAS.M$beta,
              sd_target = Exp.GWAS.F$se,
              sd_other = Exp.GWAS.M$se,
              method = "meta")

## APM
APM.F <- SEB(target_effect = Exp.GWAS.F$beta,
             other_effect = Exp.GWAS.M$beta,
             sd_target = Exp.GWAS.F$se,
             sd_other = Exp.GWAS.M$se,
             method = "APM")

## freq AW
AW.F <- AW.Est.fn(target_effect = Exp.GWAS.F$beta,
                  other_effect = Exp.GWAS.M$beta,
                  sd_target = Exp.GWAS.F$se,
                  sd_other = Exp.GWAS.M$se)


#####################################################
## compute female-specific causal effect estimate
#####################################################
## variant-exposure: raw
res.F <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                             Outcome_se = outcome.GWAS$se,
                             Exposure_effect = Exp.GWAS.F$beta,
                             Exposure_se = Exp.GWAS.F$se,
                             N = Input_cML)


## variant-exposure: meta
res.meta.F <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                                  Outcome_se = outcome.GWAS$se,
                                  Exposure_effect = meta.F$Est,
                                  Exposure_se = meta.F$sd,
                                  N = Input_cML)

## variant-exposure: APM
res.APM.F <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                                 Outcome_se = outcome.GWAS$se,
                                 Exposure_effect = APM.F$Est,
                                 Exposure_se = APM.F$sd,
                                 N = Input_cML)


## variant-exposure: AW
res.AW.F <- MR_methods_function(Outcome_effect = outcome.GWAS$beta,
                                Outcome_se = outcome.GWAS$se,
                                Exposure_effect = AW.F$AW.Est,
                                Exposure_se = AW.F$sd.AW.Est,
                                N = Input_cML)


#####################################################
## finalize results
#####################################################
Final.res <- rbind(res.F, res.meta.F, res.APM.F, res.AW.F)

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

