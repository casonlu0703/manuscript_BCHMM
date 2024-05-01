#############################################################
# MANUSCRIPT: A Bayesian Circadian Hidden Markove Model
# to Infer Rest-Activity Rhythms Using 24-hour Actigraphy Data.
# R SCRIPT ----------------------
# 1. Fit BCHMM using Stan code
#############################################################


# set up ------------------------------------------------------------------
#rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

if (!require(dplyr)) install.packages("dplyr")
if (!require(rstan)) install.packages("rstan")
# MCMC diagnostic plot
if(!require(bayesplot)) install.packages("bayesplot")
if(!require(ggplot2)) install.packages("ggplot2")

# source("../simulation/code/function_model fit_stan.R")

# fit stan model for each subject
run_stan <- function(run.id){
  
  dat <- readRDS("data/sample_data.rds") %>% 
    filter(SEQN %in% run.id) %>%
    dplyr::select(SEQN, fivemin, activity_avg)
  
  dat$y <- sqrt(dat$activity_avg)
  
  sf <- 1/12 #sampling frequency in hour unit, i.e. here sf=5min/60min=1/12
  dat$lag <- seq(1, nrow(dat))
  dat$sin_part<-sin(2*pi*dat$lag*sf/24)
  dat$cos_part<-cos(2*pi*dat$lag*sf/24)
  
  # user-defined initial values for mean and variance parameters in the Gaussian emission distributions
  state.mean <- c(2,6,10)
  state.var  <- c(2,2,2)^2
  
  # fit model using Stan
  nstate <- 3 # number of hidden states
  ncoef <- 3  # number of coefficients in the tpm
  
  stan.init <- function() {
    list(mu = state.mean, sigma = state.var
         , beta = array(
            rep(0, nstate*(nstate-1)*ncoef)
            , dim = c(nstate*(nstate-1),ncoef)
           )
         , pi1 = rep(1/3, 3))
  }           
  
  stan_data <- list(N = nrow(dat),
                    K = nstate,
                    M = ncoef, # 2 components plus a vector of 1s
                    y = dat$y, # observations
                    u = as.array(cbind(rep(1, nrow(dat)), dat$sin_part, dat$cos_part)),
                    # priors for mu and sigma
                    mu0 = mean(dat$y, na.rm = T), k0 = 1,
                    v0 = 2, sigma02 =0.5
                    #, alpha = c(1, 1, 1)
                    )
  
  fit <- stan(file = "function_stan_for_bchmm.stan" # stan code
              , data = stan_data
              , iter = 5000, chains = 1
              , init = stan.init
  )

  saveRDS(fit, paste0("results/","bchmm_stanfit_id_",run.id,".rds"))
}

# read in data ------------------------------------------------------------

# subject id
id <- unique(readRDS("data/sample_data.rds")$SEQN)

# 1.model fitting (1st fit) -------------------------------------------------

for (i in 1:length(id)){ # loop through subject IDs
  
  run.id <- id[i]
  cat(paste0("i=", i, ",ID=", run.id, sep = ""))
  
  run_stan(run.id)
}
