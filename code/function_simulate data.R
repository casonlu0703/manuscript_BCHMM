#### Simulate data ####

### simulate one data ###
simu.data <- function(mu=mu                  # a vector of mean
                      , sigma2y=sigma2y      # a vector of variance
                      , p0=prob_state_t0[1,] # a vector initial probabilities
                      , tpm_t=tpm_t          # a array of transition probabilities
                      , nstate=nstate        # number of states
                      , nobs=length(xt))     # number of observations
{   
  # generate hidden markov chain
  
  ## hidden markov chain
  mc <- rep(NA, nobs)
  ## observations
  obs <- rep(NA, nobs)
  
  ## initial state
  mc[1] <- sample(x = seq(1, nstate, by = 1), 1, prob = p0)
  for (k in 1:nstate){
    if (mc[1]==k) obs[1] <- rnorm(1, mean=mu[k], sd=sqrt(sigma2y[k]))
  }
  
  ## current state depends on the immediate one previous state only
  for (i in 2:nobs) {
    # transition probability matrix
    gamma <- tpm_t[,,i-1]
    # state at time i
    mc[i] <- sample(x = seq(1, nstate, by = 1), 1, prob = gamma[mc[i-1], ])
    # observation value at time i
    for (k in 1:nstate){
      if (mc[i]==k) obs[i] <- rnorm(1, mean=mu[k], sd=sqrt(sigma2y[k]))
    }
  }
  return(list(obs, mc))
}
