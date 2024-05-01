#####################################################################
# Simulate data and fit BCHMM using Stan
#####################################################################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(depmixS4);library(rstan);library(ggplot2);library(dplyr)


# function to simulate data
source("function_simulate data.R")

# adjust memory limit to run Stan
# memory.limit(size = 60000) 

stan.onerun <- function(data, nstate, ncoef, stan.init="random"
                        , nchains = 2, niter = 5000){
  data1 = data
  stan_data <- list(N = nrow(data1),
                    K = nstate,
                    M = ncoef, # 2 components plus a vector of 1s
                    y = data1$obs, # observations
                    u = as.array(cbind(rep(1, nrow(data1)), data1$sin_part, data1$cos_part)),
                    # priors for mu and sigma
                    mu0 = .1, k0 = .1,
                    v0 = 1, sigma02 =1
  )
  hmm_fit <- stan(file = "function_stan_for_bchmm.stan"
                  , data = stan_data
                  , iter = niter, thin = 1, chains = nchains
                  , init = stan.init)
  return(hmm_fit) # return the model fitting object
}


simu <- function(scenario, nrun){
  
  # number of hidden states
  nstate <- 3
  # number of circadian oscillations, plus the intercept
  ncoef <- 3
  
  prob_state_t <- readRDS("true_param.rds")[[1]]
  tpm_t <- readRDS("true_param.rds")[[2]]
  sin_part <- readRDS("true_param.rds")[[3]]
  cos_part <- readRDS("true_param.rds")[[4]]
  
  # scenarios, differ by mean and sd in gaussian emission dist'n
  if (scenario==1) {
    mu <- c(2, 6, 11)
    sigma2y <- c(.5, .5, .5)
  } else if (scenario==2) {
    mu <- c(2, 6, 11)
    sigma2y <- c(1, 1, 1)
  } else if (scenario==3) {
    mu <- c(2, 6, 11)
    sigma2y <- c(.6, 1.75, 2.2)
  }
  
  # simulate data
  dt <- simu.data(mu=mu                  # a vector of mean
                  , sigma2y=sigma2y      # a vector of variance
                  , p0=prob_state_t[1,]  # a vector initial probabilities
                  , tpm_t=tpm_t          # a array of transition probabilities
                  , nstate=nstate        # number of states
                  , nobs=length(sin_part)# number of observations
  )
  
  saveRDS(dt, paste0("results/simu_data_","scenario",scenario,"_",nrun,".rds"))
  
  data0 <- data.frame(obs=dt[[1]], state=dt[[2]], sin_part=sin_part, cos_part=cos_part)[-c(1:288),]

  fit2 <- stan.onerun(data0, nstate, ncoef)
  saveRDS(fit2, paste0("results/","stan_"  ,"scenario",scenario,"_",nrun,".rds"))
  rm(fit2)
  gc()
}

for(i in 1:100){
  gc()
  simu(scenario=3,nrun=i)
  Sys.sleep(100)
  gc()
}

for(i in 1:100){
  gc()
  simu(scenario=2,nrun=i)
  Sys.sleep(100)
  gc()
}

for(i in 1:100){
  gc()
  simu(scenario=3,nrun=i)
  Sys.sleep(100)
  gc()
}
