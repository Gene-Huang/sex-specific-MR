###########################################################
## semi-empirical Bayes function
###########################################################
###########################################################
## function: SEB()
## function for computing FE meta and APM estimate
###########################################################
## inputs: 
## input1: target effect size
## input2: other effect size
## input3: standard error of target effect size 
## input4: standard error of other effect size
## input5: method = c("meta", "APM")
## Notes: all these inputs can be vector format
###########################################################
## output is a "list" containing three objects
## output1: posterior mean
## output2: posterior sd
## output3: p-value of estimates
###########################################################

## method = c("meta", "APM")
SEB <- function(target_effect, other_effect, sd_target, sd_other, method = c("meta", "APM")){
  
  ###################################
  ## fixed effect meta estimation
  ###################################
  ## inverse variance of target effect group
  inv_var_target <- 1/(sd_target^2)
  
  ## inverse variance of other effect group
  inv_var_other <- 1/(sd_other^2)
  
  ## compute the weighting parameter of meta estimate
  Meta.weight <- inv_var_target/(inv_var_target + inv_var_other)
  
  ## meta estimate
  Meta.Est <- (inv_var_target*target_effect + inv_var_other*other_effect)/(inv_var_target + inv_var_other)
  
  ## variance of meta estimate
  Meta.var <- 1/(inv_var_target + inv_var_other)
  
  Meta.sd <- sqrt(Meta.var)
  
  Meta.PValue <- pchisq((Meta.Est/Meta.sd)^2, 1, lower.tail = F)
  
  if (method == "meta"){
  
  return(as.list(data.frame(Est = Meta.Est, sd = Meta.sd, pvalue = Meta.PValue)))
    
  } else {
    
    ###################################
    ## APM estimate
    ###################################
    prior_var <- ((target_effect - other_effect)^2) + (sd_other^2)

    Data_var <- sd_target^2

    ## compute the weighting parameter
    weight <- (prior_var)/(Data_var + prior_var)
    
    ## APM estimate
    posterior.mean <- (weight*target_effect) + ((1 - weight)*Meta.Est)

    ## variance of APM estimate
    posterior.variance <- (Data_var*prior_var)/(Data_var + prior_var)
    
    posterior.sd <- sqrt(posterior.variance)
    
    posterior.PValue <- pchisq((posterior.mean/posterior.sd)^2, 1, lower.tail = F)
    
    ## outputs: posterior mean, sd, p-value
    return(as.list(data.frame(Est = posterior.mean, sd = posterior.sd, pvalue = posterior.PValue)))

  }
  
}

