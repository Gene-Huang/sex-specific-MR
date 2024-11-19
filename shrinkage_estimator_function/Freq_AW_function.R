###########################################################
## compute adaptive weight (AW) estimate
###########################################################
###########################################################
## function: AW.Est.fn()
## function for computing adaptive weight estimate
###########################################################
## inputs: 
## input1: target effect size
## input2: other effect size
## input3: standard error of target effect size 
## input4: standard error of other effect size
## Notes: all these inputs can be a vector format
###########################################################
## output is a "list" containing three objects
## output1: adaptive weight estimate
## output2: standard error of adaptive weight estimate
## output3: p-value of adaptive weight estimate
###########################################################


## adaptive weight estimate function

AW.Est.fn <- function(target_effect, other_effect, sd_target, sd_other){
  
## step1: compute meta estimate (standard inverse variance weighting approach)
  
  ## inverse variance of target effect group
  inv_var_target <- 1/(sd_target^2)
  
  ## inverse variance of other effect group
  inv_var_other <- 1/(sd_other^2)
  
  ## compute the weighting parameter of meta estimate
  ## we need this value when calculating the variance of adaptive weight estimate later
  Meta.weight <- inv_var_target/(inv_var_target + inv_var_other)
  
  ## the meta estimate
  Meta.Est <- (inv_var_target*target_effect + inv_var_other*other_effect)/(inv_var_target + inv_var_other)

## step2: compute the adaptive weight estimate (AW.Est)
  
  ## compute the difference of beta effect size between two groups
  ## notation "delta" in our draft
  Effect.Diff <- target_effect - other_effect
  
  ## compute the adaptive weight parameter
  AD.weight <- (Effect.Diff^2)/((Effect.Diff^2) + (sd_target^2))
  
  ## the adaptive weight estimate
  AW.Est <- (AD.weight*target_effect) + ((1-AD.weight)*Meta.Est)
  
## step3: compute variance of AW.Est
  
  ## due to the formula is too complex
  ## I divided it into numerator and denominator part 
  ## formula can be found in Li et al., 2010 (Genetic Epidemiology)
  
  ## numerator part of the variance formula
  numerator <- (sd_target^2)*(1-Meta.weight)*(Effect.Diff^2-sd_target^2)
  
  ## denominator part of the variance formula
  denominator <- ((sd_target^2) + Effect.Diff^2)^2
  
  ## the variance of AW.Est
  Var.AW.Est <- ((sd_target^2)*((1 + numerator/denominator)^2)) + ((sd_other^2)*((numerator/denominator)^2))
  
  sd.AW.Est <- sqrt(Var.AW.Est)
  
  P.Value <- pchisq((AW.Est/sd.AW.Est)^2, 1, lower.tail = F)
  
## outputs: adaptive weight estimate and its variance
  return(as.list(data.frame(AW.Est, sd.AW.Est, P.Value)))

}




