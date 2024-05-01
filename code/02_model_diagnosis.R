#############################################################
# MANUSCRIPT: A Bayesian Circadian Hidden Markove Model
# to Infer Rest-Activity Rhythms Using 24-hour Actigraphy Data.
# R SCRIPT ----------------------
# 2. Model diagnostics
# 2.1 Continue running the model if convergence failed.
# 2.2 Check diagnostic statistics of Gelman-Rubin statistics and 
#     effective sample size.
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

# functions ---------------------------------------------------------------


# diagnostics plots
diagnosis_stan <- function(stanfit){
  
  fit2 <- stanfit
  rm(stanfit)
  
  # extract Rhat. If the largest greater than 1.05, convergence may not be good
  
  # Rhat
  p1 <- mcmc_rhat(rhat = rhat(fit2, pars = c("mu","sigma2y","beta"))) + yaxis_text(hjust = 0)
  
  # Effective size
  ratios_cp <- neff_ratio(fit2, pars = c("mu","sigma2y","beta"))
  p2 <- mcmc_neff(ratios_cp, size = 2) + yaxis_text(hjust = 0)
  
  #extract a iterations x chains x parameters array of posterior draws
  posterior_cp <- as.array(fit2)
  # log of posterior density
  lp_cp <- log_posterior(fit2)
  # diagnostic parameters from NUTS
  np_cp <- nuts_params(fit2)
  
  # Trace plot
  p3 <- mcmc_trace(posterior_cp, regex_pars = "^mu\\[", np=np_cp) +  xlab("Post-warmup iteration")
  p4 <- mcmc_trace(posterior_cp, regex_pars = "^sigma2y\\[", np=np_cp) + xlab("Post-warmup iteration")
  p5 <- mcmc_trace(posterior_cp, regex_pars = "^beta\\[", np=np_cp) + xlab("Post-warmup iteration")
  p6 <- mcmc_nuts_divergence(np_cp, lp_cp, chain = 1)
  #p7 <- mcmc_nuts_divergence(np_cp, lp_cp, chain = 2)
  
  # posterior distribution
  p8 <- mcmc_areas(as.matrix(fit2), regex_pars = "mu", prob = 0.8)
  p9 <- mcmc_areas(as.matrix(fit2), regex_pars = "sigma2y", prob = 0.8)
  p10 <- mcmc_areas(as.matrix(fit2), regex_pars = "beta", prob = 0.8)
  
  print(p1);print(p2);print(p3);print(p4);print(p5)
  print(p6); #print(p7);
  print(p8);print(p9);print(p10)
  
}


# perform bchmm for subjects with not converged stanfit in the 1st run
# use Depmix results as initial values
rerun_stan <- function(run.id){
  
  dat <- readRDS("data/sample_data.rds") %>% 
    filter(SEQN %in% run.id) %>% 
    dplyr::select(SEQN, fivemin, activity_avg)
  
  dat$y <- sqrt(dat$activity_avg)
  
  sf <- 1/12 #sampling frequency in hour unit, i.e. here sf=5min/60min=1/12
  dat$lag <- seq(1, nrow(dat))
  dat$sin_part<-sin(2*pi*dat$lag*sf/24)
  dat$cos_part<-cos(2*pi*dat$lag*sf/24)
  
  #### fit stan #####
  nstate <- 3 # number of hidden states
  ncoef <- 3  # number of coefficients in the tpm
  
  stan_data <- list(N = nrow(dat),
                    K = nstate,
                    M = ncoef, # 2 components plus a vector of 1s
                    y = dat$y, # observations
                    u = as.array(cbind(rep(1, nrow(dat)), dat$sin_part, dat$cos_part)),
                    # priors for mu and sigma
                    mu0 = mean(dat$y, na.rm = T), k0 = 1,
                    v0 = 2, sigma02 =0.5)
  
  fit0 <- readRDS(paste0("results/stanfits/","bchmm_stanfit_id_",run.id,".rds"))
  
  # discard the not converged current model fit
  saveRDS(fit0, paste0("results/stanfits/discard/","bchmm_stanfit_id_",run.id,"_notconvg.rds"))
  
  # number of iterations
  niter <- 10000
  
  # identify convergence by Rhat. 
  # if max Rhat > 1.05, OR if min effective size ratio < .5
  # rerun the stan model fit with better initial values
  max.Rhat <- max(rhat(fit0))
  min.neff <- min(neff_ratio(fit0, pars = c("mu","sigma2y","beta")))*5000/2
  
  # posterior means from last fit as initial values
  stan.init <- get_posterior_mean(fit0)
  
  while (max.Rhat > 1.05 | min.neff < 100 ) {
    
    fit <- stan(file = "function_stan_for_bchmm.stan"
                , data = stan_data
                , iter = niter, chains = 1
                , init = stan.init
    )
    max.Rhat <- max(rhat(fit))
    min.neff <- min(neff_ratio(fit, pars = c("mu","sigma2y","beta")))*niter
    
    rm(stan.init)
    stan.init <- get_posterior_mean(fit)
    
    # replace with the new stanfit
    rm(fit0); fit0 <- fit
    rm(fit)
    
  }
  
  # once model converged, remove previous fit from current folder
  file.remove(paste0("results/stanfits/","bchmm_stanfit_id_",run.id,".rds"))
  # save the converged model fit
  saveRDS(fit0, paste0("results/stanfits/","bchmm_stanfit_id_",run.id,".rds"))
}

# read in data ------------------------------------------------------------

# subject id
id <- unique(readRDS("data/sample_data.rds")$SEQN)

# 2.1 model diagnostics -------------------------------------------------------

rerun.id <- c()
j <- 1

for (i in 1:length(id)){
  run.id <- id[i]
  
  #pdf(paste0("results/bchmm_diagnosis_id_",run.id,"_rerun.pdf"))
  fit2 <- readRDS(paste0("results/stanfits/","bchmm_stanfit_id_",run.id,".rds")) 
  #diagnosis_stan(stanfit = fit2)
  #dev.off()  
  
  # identify convergence by Rhat. 
  # if max Rhat > 1.05, OR if min effective size ratio < .5
  # rerun the stan model fit with better initial values
  max.Rhat <- max(rhat(fit2))
  min.neff <- min(neff_ratio(fit2, pars = c("mu","sigma2y","beta")))
  
  #if (max.Rhat>1.05 | min.neff < 0.5 ){ # if not converge, rerun
  if (min.neff*2500 <= 100 ){ # if not converge, rerun
    print(paste0("Reruning! Stanfit for ", i,"-th subject, ID=", run.id, " is not converge."))
    rerun.id[j] <- run.id; j <- j+1
    
  } else { # if converge, skip
    print(paste0("Congrats! Stanfit for ", i,"-th subject, ID=", run.id, " converged."))
  }
  rm(list = c("fit2", "max.Rhat", "min.neff"))
}

#saveRDS(rerun.id, "rerun_id_neff_lt_0.1.rds")
saveRDS(rerun.id, "rerun_id_neff_lt_100.rds")

# 2.2. Rerun not converged model fit -------------------------------------------

rerun.id <- readRDS("rerun_id_neff_lt_100.rds")

for (i in 1:length(rerun.id)){
  run.id <- rerun.id[i]
  cat(paste0("i=", i, ",ID=", run.id, sep = ""))
  rerun_stan(run.id)
}
