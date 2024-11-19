###############################################################
## function for applying two-sample MR methods
###############################################################
###############################################################
## function: MR_methods_function()
###############################################################
## selected methods: 
## weighted median
## IVW, penalized and robust IVW
## MR-Egger, penalized and robust MR-Egger
## MR-PRESSO
## MR Contamination mixture 
## Constrained maximum likelihood (cML) 
## MR-RAPS (with robust Tukey's loss)
## Note: if the number of IVs is less than or equal to 3
## we just output the result from IVW method
###############################################################
###############################################################
## inputs: 
## input1: variant-outcome association effect size
## input2: standard error of variant-outcome effect size 
## input3: variant-exposure association effect size
## input4: standard error of variant-exposure effect size
## input5: sample size of GWAS (the smallest sample size between exposure and outcome is recommended)
## Notes: input 5 is for applying "cML" method
###############################################################
###############################################################
## load R packages first
###############################################################
library(dplyr)

## load MR-PRESSO R package
library(MRPRESSO)

## load mr.raps R package
library(mr.raps)

## load TwoSampleMR R package
library(MendelianRandomization)
###############################################################


MR_methods_function <- function(Outcome_effect, Outcome_se, Exposure_effect, Exposure_se, N){
  
  ## if the input does not contain enough IVs (less than or equal to 3)
  ## just output the result from IVW
  if(length(Exposure_effect) < 4){
    MR_results <- mr_allmethods(mr_input(bx = Exposure_effect,
                                         bxse = Exposure_se,
                                         by = Outcome_effect,
                                         byse = Outcome_se),
                                method = "ivw", iterations = 10000)
    ## select IVW result
    Index <- which(MR_results$Values[,"Method"] == "IVW")
    
    Final_res <- MR_results$Values[Index, ]
    
    colnames(Final_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue")
    
  } else{
    
    ############################
    ## MR-PRESSO
    ############################
    Data <- data.frame(Y_effect = Outcome_effect, 
                       Y_se = Outcome_se,
                       E_effect = Exposure_effect,
                       E_se = Exposure_se)
    
    MR_PRESSO <- mr_presso(BetaOutcome = "Y_effect",
                           BetaExposure = "E_effect",
                           SdOutcome = "Y_se",
                           SdExposure = "E_se",
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                           data = Data, NbDistribution = 1000,  SignifThreshold = 0.05)
    
    ## if no outlier detect
    ## MR-PRESSO outputs IVW estimation
    if (is.na(MR_PRESSO$`Main MR results`[2,]$`Causal Estimate`) == T){
      MR_PRESSO_res <- select(MR_PRESSO$`Main MR results`[1,], "Causal Estimate", "Sd", "P-value")
    } else {
      MR_PRESSO_res <- select(MR_PRESSO$`Main MR results`[2,], "Causal Estimate", "Sd", "P-value")
    }
    
    MR_PRESSO_res <- data.frame(Method = "MR_PRESSO", 
                                Estimate = MR_PRESSO_res$`Causal Estimate`, 
                                se = MR_PRESSO_res$Sd,
                                CI_L = MR_PRESSO_res$`Causal Estimate` - 1.96*MR_PRESSO_res$Sd,
                                CI_U = MR_PRESSO_res$`Causal Estimate` + 1.96*MR_PRESSO_res$Sd,
                                Pvalue = MR_PRESSO_res$`P-value`)

    ############################
    ## MR-RAPS with robust Tukey's loss
    ## use overdispersion version
    ############################
    MR_RAPS_robust <- mr.raps.overdispersed.robust(b_exp = Exposure_effect,
                                                   b_out = Outcome_effect,
                                                   se_exp = Exposure_se,
                                                   se_out = Outcome_se,
                                                   loss.function = "tukey", 
                                                   niter = 10000,
                                                   suppress.warning = TRUE)
    
    
    MR_RAPS_robust_res <- data.frame(Method = "MR_RAPS",
                                     Estimate = MR_RAPS_robust$beta.hat, 
                                     se = MR_RAPS_robust$beta.se,
                                     CI_L = MR_RAPS_robust$beta.hat - 1.96*MR_RAPS_robust$beta.se,
                                     CI_U = MR_RAPS_robust$beta.hat + 1.96*MR_RAPS_robust$beta.se,
                                     Pvalue = MR_RAPS_robust$beta.p.value)
    
    ############################
    ## other MR methods
    ## IVW, robust IVW, MR-Egger, robust MR-Egger
    ############################
    MR_main <- mr_allmethods(mr_input(bx = Exposure_effect,
                                      bxse = Exposure_se,
                                      by = Outcome_effect,
                                      byse = Outcome_se),
                             method = "all", iterations = 10000)
    
    ## select results from: weighted median, IVW, robust IVW, ME-Egger, robust ME-Egger
    Index <- which(MR_main$Values[,"Method"] %in% c("Weighted median", "IVW", "Penalized robust IVW", "MR-Egger", "Penalized robust MR-Egger"))
    MR_main_res <- MR_main$Values[Index, ]
    colnames(MR_main_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue")
    
    
    ############################
    ## MR Contamination mixture 
    ############################
    MR_mix <- mr_conmix(mr_input(bx = Exposure_effect,
                                 bxse = Exposure_se,
                                 by = Outcome_effect,
                                 byse = Outcome_se),
                        psi = 0, CIMin = -3, CIMax = 3, CIStep = 0.001)
    
    ## note: MR-mixture cannot output se
    MR_mix_res <- data.frame(MR_mix$class[1],
                             MR_mix$Estimate,
                             NA,
                             MR_mix$CILower,
                             MR_mix$CIUpper,
                             MR_mix$Pvalue)
    
    MR_mix_res <- MR_mix_res[1, ]
    
    colnames(MR_mix_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue")
    
    
    ############################
    ## MR Constrained maximum likelihood 
    ############################
    MR_CML <- mr_cML(mr_input(bx = Exposure_effect,
                              bxse = Exposure_se,
                              by = Outcome_effect,
                              byse = Outcome_se), 
                     maxit = 100, MA = TRUE, DP = FALSE, n = N)
    
    MR_CML_res <- data.frame(MR_CML$class[1],
                             MR_CML$Estimate,
                             MR_CML$StdError,
                             MR_CML$CILower,
                             MR_CML$CIUpper,
                             MR_CML$Pvalue)
    
    colnames(MR_CML_res) <- c("Method", "Estimate", "se", "CI_L", "CI_U", "Pvalue")
    
    ############################
    ## combine all results
    ############################
    Final_res <- rbind(MR_main_res,
                       MR_PRESSO_res,
                       MR_mix_res,
                       MR_CML_res, 
                       MR_RAPS_robust_res)
    
  }
    
  return(Final_res)
  
}






