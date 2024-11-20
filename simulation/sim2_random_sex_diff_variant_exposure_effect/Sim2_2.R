##############################################################
## Simulation setting 2-2
## generate IV effect by bi-variate normal distribution
## mean of exposure GWAS effect size = 0.1 (both female and male)
## fix correlation between female and male G to E effect size as 0.7
## different causal effect for female and male (female = 0, 0.05, 0.15; male = 0.1)
## competing methods:
## weighted median, IVW, penalized and robust IVW
## MR-Egger, penalized and robust MR-Egger
## mixture, cML, MR-RAPS
##############################################################
## load data table R package
## use to summarize simulation results
library(data.table)

library(MASS)

## load MR function
source("TwoSampleMR_methods_function.R")

## load semi-empirical Bayes function
source("SEB_function.R")

## load freq AW function
source("Freq_AW_function.R")

## number of replications
RP <- 1000

## sample size of population 1
## generate exposure GWAS
P1.Sample.F <- 2000

P1.Sample.M <- 20000


## sample size of population 2
## generate outcome GWAS
P2.Sample.F <- 10000

P2.Sample.M <- 10000


## set causal effect: exposure to outcome
## sim2-2: different causal effect
## causal effect female = 0, 0.05, 0.15
## causal effect male = 0.1
True.Causal.F <- 0
#True.Causal.F <- 0.05
#True.Causal.F <- 0.15

True.Causal.M <- 0.1

True.Causal <- c(Female = True.Causal.F, Male = True.Causal.M)

## number of SNPs (IV)
N.IV <- 100

## allele freq
## set the same for all SNPs
af <- 0.3

## construct G to E effect for both female and male
## N.IV*2 matrix
## first column: female G to E effect size
## second column: male G to E effect size
Sim.G.to.E <- matrix(0, N.IV, 2)
colnames(Sim.G.to.E) <- c("Female", "Male")

## consider the correlation of G to E effect sizes between female and male (across all IVs)
Corr <- 0.7

VAR <- 0.01

## covariance
COV <- Corr*sqrt(VAR)*sqrt(VAR)

## the mean value of G to E effect size
mean.G.to.E <- c(0.1, 0.1)

## the covariance matrix to generate G to E effect size
## a 2*2 matrix
Cov.matrix <- matrix(c(VAR, COV, COV, VAR), 2, 2)

## set the true effect size from unknown confounding to exposure
## assume the same across female and male
Sim.U.to.E <- 0.1

## set the true effect size from unknown confounding to outcome
## assume the same across female and male
Sim.U.to.O <- 0.1


## save the marginal summary statistics (sex-specific exposure/outcome GWAS results)
F.Outcome.GWAS <- matrix(0, N.IV, 3)
colnames(F.Outcome.GWAS) <- c("Beta", "se", "pvalue")

F.Exposure.GWAS <- matrix(0, N.IV, 3)
colnames(F.Exposure.GWAS) <- c("Beta", "se", "pvalue")

M.Outcome.GWAS <- matrix(0, N.IV, 3)
colnames(M.Outcome.GWAS) <- c("Beta", "se", "pvalue")

M.Exposure.GWAS <- matrix(0, N.IV, 3)
colnames(M.Exposure.GWAS) <- c("Beta", "se", "pvalue")


####################################################
## save the causal estimate from each competing methods
####################################################
## for female
####################################################
## raw female
res.F <- vector(mode = "list", length = RP)

## meta estimate
res.meta.F <- vector(mode = "list", length = RP)

## APM
res.APM.F <- vector(mode = "list", length = RP)

## freq AW
res.AW.F <- vector(mode = "list", length = RP)


####################################################
## for male
####################################################
## raw male
res.M <- vector(mode = "list", length = RP)

## meta estimate
res.meta.M <- vector(mode = "list", length = RP)

## APM
res.APM.M <- vector(mode = "list", length = RP)

## freq AW
res.AW.M <- vector(mode = "list", length = RP)


####################################################
## start simulation studies
####################################################

for (i in 1:RP){
  
  ## generate female/male G to E effect
  ## from bi-variate normal distribution
  ## each variant was iid sampled
  IV.effect <- mvrnorm(n = N.IV, mu = mean.G.to.E, Sigma = Cov.matrix)
  
  Sim.G.to.E[, "Female"] <- IV.effect[,1]
  
  Sim.G.to.E[, "Male"] <- IV.effect[,2]
  
  ## population 1
  ## generate SNP data
  P1.G.F <- matrix(rbinom(P1.Sample.F*N.IV, 2, af), P1.Sample.F, N.IV)
  
  P1.G.M <- matrix(rbinom(P1.Sample.M*N.IV, 2, af), P1.Sample.M, N.IV)
  
  ## generate unknown confounding
  ## for U to E
  ## female and male have "different" distribution (Normal dist. with same mean but different variance)
  P1.U.F <- rnorm(P1.Sample.F, 0, 1)
  
  P1.U.M <- rnorm(P1.Sample.M, 0, 0.5)
  
  ## generate exposure 
  P1.Exposure.F <- 1 + Sim.G.to.E[,"Female"] %*% t(P1.G.F) + Sim.U.to.E*P1.U.F + rnorm(P1.Sample.F, 0, 1)
  
  P1.Exposure.M <- 1 + Sim.G.to.E[,"Male"] %*% t(P1.G.M) + Sim.U.to.E*P1.U.M + rnorm(P1.Sample.M, 0, 1)
  
  
  ## population 2
  ## generate SNP data
  P2.G.F <- matrix(rbinom(P2.Sample.F*N.IV, 2, af), P2.Sample.F, N.IV)
  
  P2.G.M <- matrix(rbinom(P2.Sample.M*N.IV, 2, af), P2.Sample.M, N.IV)
  
  ## generate unknown confounding
  ## for U to E
  ## female and male have "different" distribution
  P2.U.F <- rnorm(P2.Sample.F, 0, 1)
  
  P2.U.M <- rnorm(P2.Sample.M, 0, 0.5)
  
  ## generate exposure 
  ## note: the G to E effect size set the same for population 1 and 2
  P2.Exposure.F <- 1 + Sim.G.to.E[,"Female"] %*% t(P2.G.F) + Sim.U.to.E*P2.U.F + rnorm(P2.Sample.F, 0, 1)
  
  P2.Exposure.M <- 1 + Sim.G.to.E[,"Male"] %*% t(P2.G.M) + Sim.U.to.E*P2.U.M + rnorm(P2.Sample.M, 0, 1)
  
  ## generate outcome
  Outcome.F <- 1 + True.Causal["Female"]*P2.Exposure.F + Sim.U.to.O*P2.U.F + rnorm(P2.Sample.F, 0, 1)
  
  Outcome.M <- 1 + True.Causal["Male"]*P2.Exposure.M + Sim.U.to.O*P2.U.M + rnorm(P2.Sample.M, 0, 1)
  
  ## prepare the data format for running linear regression
  Ex.F <- cbind(Exposure.F = as.numeric(P1.Exposure.F), P1.G.F)
  Ex.F <- as.data.frame(Ex.F)
  
  Out.F <- cbind(Outcome.F = as.numeric(Outcome.F), P2.G.F)
  Out.F <- as.data.frame(Out.F)
  
  Ex.M <- cbind(Exposure.M = as.numeric(P1.Exposure.M), P1.G.M)
  Ex.M <- as.data.frame(Ex.M)
  
  Out.M <- cbind(Outcome.M = as.numeric(Outcome.M), P2.G.M)
  Out.M <- as.data.frame(Out.M)
  
  
  ## fit model
  ## get marginal summary statistics
  ## Exposure GWAS
  ## Outcome GWAS
  for (j in 1:N.IV){
    
    ## female outcome GWAS
    F.Outcome.GWAS[j, ] <- coef(summary(lm(Outcome.F ~ Out.F[, (j + 1)], data = Out.F)))[2, c("Estimate", "Std. Error", "Pr(>|t|)")]
    
    ## female exposure GWAS
    F.Exposure.GWAS[j, ] <- coef(summary(lm(Exposure.F ~ Ex.F[, (j + 1)], data = Ex.F)))[2, c("Estimate", "Std. Error", "Pr(>|t|)")]
    
    ## male outcome GWAS
    M.Outcome.GWAS[j, ] <- coef(summary(lm(Outcome.M ~ Out.M[, (j + 1)], data = Out.M)))[2, c("Estimate", "Std. Error", "Pr(>|t|)")]
    
    ## male exposure GWAS
    M.Exposure.GWAS[j, ] <- coef(summary(lm(Exposure.M ~ Ex.M[, (j + 1)], data = Ex.M)))[2, c("Estimate", "Std. Error", "Pr(>|t|)")]
    
  }
  
  
  ###########################
  ## adjusted variant-exposure effect estimate for female
  ###########################
  meta.F <- SEB(F.Exposure.GWAS[,"Beta"], M.Exposure.GWAS[,"Beta"], F.Exposure.GWAS[,"se"], M.Exposure.GWAS[,"se"], method = "meta")
  
  APM.F <- SEB(F.Exposure.GWAS[,"Beta"], M.Exposure.GWAS[,"Beta"], F.Exposure.GWAS[,"se"], M.Exposure.GWAS[,"se"], method = "APM")
  
  AW.F <- AW.Est.fn(F.Exposure.GWAS[,"Beta"], M.Exposure.GWAS[,"Beta"], F.Exposure.GWAS[,"se"], M.Exposure.GWAS[,"se"])
  
  ###########################
  ## adjusted variant-exposure effect estimate for male
  ###########################
  meta.M <- SEB(M.Exposure.GWAS[,"Beta"], F.Exposure.GWAS[,"Beta"], M.Exposure.GWAS[,"se"], F.Exposure.GWAS[,"se"], method = "meta")
  
  APM.M <- SEB(M.Exposure.GWAS[,"Beta"], F.Exposure.GWAS[,"Beta"], M.Exposure.GWAS[,"se"], F.Exposure.GWAS[,"se"], method = "APM")
  
  AW.M <- AW.Est.fn(M.Exposure.GWAS[,"Beta"], F.Exposure.GWAS[,"Beta"], M.Exposure.GWAS[,"se"], F.Exposure.GWAS[,"se"])
  
  
  ###############################
  ## start MR analysis
  ## estimate causal effect
  ###############################
  ## exposure GWAS: raw female
  ###############################
  res.F[[i]] <- MR_methods_function(Outcome_effect = F.Outcome.GWAS[,"Beta"],
                                    Outcome_se = F.Outcome.GWAS[,"se"],
                                    Exposure_effect = F.Exposure.GWAS[,"Beta"],
                                    Exposure_se = F.Exposure.GWAS[,"se"],
                                    N = P1.Sample.F)
  ###############################
  ## exposure GWAS: meta female
  ###############################
  res.meta.F[[i]] <- MR_methods_function(Outcome_effect = F.Outcome.GWAS[,"Beta"],
                                         Outcome_se = F.Outcome.GWAS[,"se"],
                                         Exposure_effect = meta.F$Est,
                                         Exposure_se = meta.F$sd,
                                         N = P1.Sample.F)
  ###############################
  ## exposure GWAS: Bayesian AW female
  ###############################
  res.APM.F[[i]] <- MR_methods_function(Outcome_effect = F.Outcome.GWAS[,"Beta"],
                                        Outcome_se = F.Outcome.GWAS[,"se"],
                                        Exposure_effect = APM.F$Est,
                                        Exposure_se = APM.F$sd,
                                        N = P1.Sample.F)
  ###############################
  ## exposure GWAS: freq AW female
  ###############################
  res.AW.F[[i]] <- MR_methods_function(Outcome_effect = F.Outcome.GWAS[,"Beta"],
                                        Outcome_se = F.Outcome.GWAS[,"se"],
                                        Exposure_effect = AW.F$AW.Est,
                                        Exposure_se = AW.F$sd.AW.Est,
                                        N = P1.Sample.F)
  
  ###############################
  ## exposure GWAS: raw male
  ###############################
  res.M[[i]] <- MR_methods_function(Outcome_effect = M.Outcome.GWAS[,"Beta"],
                                    Outcome_se = M.Outcome.GWAS[,"se"],
                                    Exposure_effect = M.Exposure.GWAS[,"Beta"],
                                    Exposure_se = M.Exposure.GWAS[,"se"],
                                    N = P1.Sample.M)
  
  ###############################
  ## exposure GWAS: meta male
  ###############################
  res.meta.M[[i]] <- MR_methods_function(Outcome_effect = M.Outcome.GWAS[,"Beta"],
                                         Outcome_se = M.Outcome.GWAS[,"se"],
                                         Exposure_effect = meta.M$Est,
                                         Exposure_se = meta.M$sd,
                                         N = P1.Sample.M)
  
  ###############################
  ## exposure GWAS: Bayesian AW male
  ###############################
  res.APM.M[[i]] <- MR_methods_function(Outcome_effect = M.Outcome.GWAS[,"Beta"],
                                        Outcome_se = M.Outcome.GWAS[,"se"],
                                        Exposure_effect = APM.M$Est,
                                        Exposure_se = APM.M$sd,
                                        N = P1.Sample.M)
  
  ###############################
  ## exposure GWAS: freq AW male
  ###############################
  res.AW.M[[i]] <- MR_methods_function(Outcome_effect = M.Outcome.GWAS[,"Beta"],
                                        Outcome_se = M.Outcome.GWAS[,"se"],
                                        Exposure_effect = AW.M$AW.Est,
                                        Exposure_se = AW.M$sd.AW.Est,
                                        N = P1.Sample.M)
  
  print(i)
  
}

###############################
## combine list results
###############################

## for female results
raw.Female <- do.call(rbind, res.F)

meta.Female <- do.call(rbind, res.meta.F)

APM.Female <- do.call(rbind, res.APM.F)

AW.Female <- do.call(rbind, res.AW.F)

## for male results
raw.Male <- do.call(rbind, res.M)

meta.Male <- do.call(rbind, res.meta.M)

APM.Male <- do.call(rbind, res.APM.M)

AW.Male <- do.call(rbind, res.AW.M)


########################################## 
## summarize results across replications
## 95% CI coverage rate
## MSE
## sd of MSE
## length of 95% CI
########################################## 
## apply data.table format
########################################## 

## female results
raw.Female <- data.table(raw.Female)
meta.Female <- data.table(meta.Female)
APM.Female <- data.table(APM.Female)
AW.Female <- data.table(AW.Female)

## male results
raw.Male <- data.table(raw.Male)
meta.Male <- data.table(meta.Male)
APM.Male <- data.table(APM.Male)
AW.Male <- data.table(AW.Male)

########################################## 
## results summary
########################################## 

## raw female results
sum.raw.Female <- raw.Female[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Female"] & True.Causal["Female"] <= CI_U, 1, 0), na.rm = T),
                                 MSE = mean((as.numeric(Estimate) - True.Causal["Female"])^2, na.rm = T),
                                 se_SE = sd(((as.numeric(Estimate) - True.Causal["Female"])^2), na.rm = T)/sqrt(RP),
                                 CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                             by = "Method"]


## meta female 
sum.meta.Female <- meta.Female[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Female"] & True.Causal["Female"] <= CI_U, 1, 0), na.rm = T),
                                   MSE = mean((as.numeric(Estimate) - True.Causal["Female"])^2, na.rm = T),
                                   se_SE = sd(((as.numeric(Estimate) - True.Causal["Female"])^2), na.rm = T)/sqrt(RP),
                                   CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                               by = "Method"]

## APM female
sum.APM.Female <- APM.Female[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Female"] & True.Causal["Female"] <= CI_U, 1, 0), na.rm = T),
                                 MSE = mean((as.numeric(Estimate) - True.Causal["Female"])^2, na.rm = T),
                                 se_SE = sd(((as.numeric(Estimate) - True.Causal["Female"])^2), na.rm = T)/sqrt(RP),
                                 CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                             by = "Method"]

## AW female 
sum.AW.Female <- AW.Female[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Female"] & True.Causal["Female"] <= CI_U, 1, 0), na.rm = T),
                                 MSE = mean((as.numeric(Estimate) - True.Causal["Female"])^2, na.rm = T),
                                 se_SE = sd(((as.numeric(Estimate) - True.Causal["Female"])^2), na.rm = T)/sqrt(RP),
                                 CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                             by = "Method"]


## raw.Male
sum.raw.Male <- raw.Male[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Male"] & True.Causal["Male"] <= CI_U, 1, 0), na.rm = T),
                             MSE = mean((as.numeric(Estimate) - True.Causal["Male"])^2, na.rm = T),
                             se_SE = sd(((as.numeric(Estimate) - True.Causal["Male"])^2), na.rm = T)/sqrt(RP),
                             CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                         by = "Method"]




## meta male 
sum.meta.Male <- meta.Male[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Male"] & True.Causal["Male"] <= CI_U, 1, 0), na.rm = T),
                               MSE = mean((as.numeric(Estimate) - True.Causal["Male"])^2, na.rm = T),
                               se_SE = sd(((as.numeric(Estimate) - True.Causal["Male"])^2), na.rm = T)/sqrt(RP),
                               CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                           by = "Method"]


## APM male 
sum.APM.Male <- APM.Male[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Male"] & True.Causal["Male"] <= CI_U, 1, 0), na.rm = T),
                             MSE = mean((as.numeric(Estimate) - True.Causal["Male"])^2, na.rm = T),
                             se_SE = sd(((as.numeric(Estimate) - True.Causal["Male"])^2), na.rm = T)/sqrt(RP),
                             CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                         by = "Method"]


## AW male
sum.AW.Male <- AW.Male[, .(Covereage = mean(ifelse(CI_L <= True.Causal["Male"] & True.Causal["Male"] <= CI_U, 1, 0), na.rm = T),
                             MSE = mean((as.numeric(Estimate) - True.Causal["Male"])^2, na.rm = T),
                             se_SE = sd(((as.numeric(Estimate) - True.Causal["Male"])^2), na.rm = T)/sqrt(RP),
                             CI.Len = mean(as.numeric(CI_U)-as.numeric(CI_L), na.rm = T)),
                         by = "Method"]


