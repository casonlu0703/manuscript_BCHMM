#############################################################
# MANUSCRIPT: A Bayesian Circadian Hidden Markove Model
# to Infer Rest-Activity Rhythms Using 24-hour Actigraphy Data.
# R SCRIPT ----------------------
# 3. calculate circadian parameters of RI, RA and CRC
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

# extract estimated parameters from stanfit
stan.param <- function(hmm_fit, nstate=3, ncoef=3){
  
  # vector of betas
  # posterior medians of circadian coefficients
  beta0.fit <- rstan::summary(hmm_fit, pars=c("beta"))$summary[,"50%"]
  
  # store into an array of matrix of coefficients corresponding to a state
  beta1.fit <- array(beta0.fit, dim=c(ncoef,nstate-1,nstate)) # dim=c(nrow, ncol, number of matrix)
  beta2.fit <- array(NA, dim=c(ncoef,nstate,nstate))
  
  for (i in 1:nstate){
    beta2.fit[,,i] <- cbind(beta1.fit[,which(c(1:(nstate-1))<=(i-1)),i]
                            , rep(0,ncoef)
                            , beta1.fit[,which(c(1:(nstate-1))>(i-1)),i])
  }

  sf <- 1/12 #sampling frequency in hour unit, i.e. here sf=5min/60min=1/12
  lag <- seq(1, 288*2)
  sin_part<-sin(2*pi*lag*sf/24)
  cos_part<-cos(2*pi*lag*sf/24)
  
  XX <- cbind(1, sin_part, cos_part)
  prob_state_t.fit <- matrix(NA, nrow = nrow(XX), ncol = nstate)
  
  # initial probabilities
  p0.fit <- get_posterior_mean(hmm_fit, pars = c("pi1"))
  prob_state_t.fit[1,] <- p0.fit
  
  # transition probabilities, last time point has no tpm
  tpm_t.fit <- array(NA, dim=c(nstate,nstate, nrow(XX)-1))
  
  for (i in 2:nrow(XX)) {
    tpm.fit <- matrix(NA, nstate, nstate)
    for (state in 1:nstate) { # time-dependent transition probability matrix
      tpm.fit[state,] <- exp(XX[i,]%*%beta2.fit[,,state])/sum(exp(XX[i,]%*%beta2.fit[,,state]))
    }
    prob_state_t.fit[i,] <- prob_state_t.fit[i-1,]%*%tpm.fit
    tpm_t.fit[,,i-1] <- tpm.fit
  }
  
  # from noon 12pm to noon 12pm
  clock.time <- seq(0,24*2-1/12,1/12)
  sp_fit2 <- prob_state_t.fit[(which(clock.time == 12)):(which(clock.time == 24+12)-1), ]
  
  return(
    list(
      param_fit2 = c(
        summary(hmm_fit, pars=c("mu"))$summary[,"50%"]
        , summary(hmm_fit, pars=c("sigma2y"))$summary[,"50%"] # variance (posterior median)
        )
      , sp_fit2 = sp_fit2 #state probabilities
      , ci = summary(hmm_fit, pars = c("pi1","mu", "sigma2y", "beta"))$summary[,c("50%", "se_mean", "sd", "2.5%", "97.5%")]
      )
    )
  
}

# read in data ------------------------------------------------------------

# subject id
id <- unique(readRDS("data/sample_data.rds")$SEQN)

# 3.extract model fits ------------------------------------------------------

files <- paste0("results/bchmm_stanfit_id_", id, ".rds")
temp <- lapply(files, readRDS)

results <- lapply(temp, stan.param)

saveRDS(results, "BCHMM_results_match_sample.rds")

results <- readRDS("BCHMM_results_match_sample.rds")

# state means and variances
state_param <- t(as.data.frame(
  lapply(results, function(x){
    return(x[[1]])
  })
))

rownames(state_param) <- NULL
colnames(state_param) <- c(paste0("state_", 1:3, "_mean"), paste0("state_", 1:3, "_var"))

# extract state probabilities
## state 1
sp1 <- t(
  as.data.frame(
    lapply(results, function(x){
      return(x[[2]][,1])
      })
    )
  )

rownames(sp1) <- NULL
colnames(sp1) <- paste0("state1_prob_", 1:288)

sp1 <- cbind(id, sp1)
  
## state 2
sp2 <- t(
  as.data.frame(
    lapply(results, function(x){
      return(x[[2]][,2])
    })
  )
)

rownames(sp2) <- NULL
colnames(sp2) <- paste0("state2_prob_", 1:288)

sp2 <- cbind(id, sp2)

## state 3
sp3 <- t(
  as.data.frame(
    lapply(results, function(x){
      return(x[[2]][,3])
    })
  )
)
rownames(sp3) <- NULL
colnames(sp3) <- paste0("state3_prob_", 1:288)
sp3 <- cbind(id, sp3)

saveRDS(list(sp1, sp2, sp3), "results/BCHMM_state_probabilities.rds")

# calculate circadian parameters
source("../simulation/code/function_calculate_circadian_parameters.R")

circadian_params <- lapply(results, function(x){
  circadian.param(one.day.prob = x[[2]])
})

RI <- unlist(lapply(circadian_params, function(x){
  x$RI
}))

RA <- unlist(lapply(circadian_params, function(x){
  x$RA
}))

CRC <- unlist(lapply(circadian_params, function(x){
  x$CRC
}))

dat.circadian <- data.frame(SEQN = id, state_param, RI = RI, RA = RA, CRC = CRC)
saveRDS(dat.circadian, "results/sample_circadian_params.rds")


# plot rest-activity profile
source("function_figure_day_profile.R")

pdf("results/rest_activity_figure.pdf")

for (i in 1:length(id)) {
  print(figure_day_profile(id = id[i]
                           , one_day_prob = results[[i]][[2]]
                           , index = seq(0,24-1/12,1/12)
                           , index_position = seq(0,24,2)
                           , index_lable = c('12','14','16','18','20','22', '0/24','2','4','6','8','10','12')
                          ))

}
dev.off()


