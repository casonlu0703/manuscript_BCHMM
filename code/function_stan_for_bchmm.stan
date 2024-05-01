data {
  int<lower=0> N; // num of observations
  int<lower=1> K; // num of states
  int<lower=0> M; // num of covariates, a column vector of 1s should be added to the covariates matrix
  
  // observed activity counts
  real y[N];
  // circadian components (or other covariates)
  matrix[N, M] u;
  // prior for mu
  real <lower=0> mu0;
  real <lower=0> k0;
  // prior for sigma
  real <lower=0> v0;
  real <lower=0> sigma02;
  // prior for delta, initial probabilities
  // vector<lower=0>[K] alpha;   // transit prior
}

parameters {
  // real mu[K];
  positive_ordered[K] mu;    // observation means
  real <lower=0> sigma2y[K]; // observation standard deviations
  
  matrix[K*(K-1), M] beta;   // coefficients for circadian components, each vector corresponding to state i transiting to state j
  simplex[K] pi1;            // initial state probabilities
  
}

model {
  
    // At each time point, there's a transition probability matrix. Each entry in a K by K matrix
    // is a linear combination of covariates (circadian components)
    matrix[K, K] unA[N]; // exponentiated linear combination in transition probabilities
    matrix[K, K] A[N];   // normalize, row sum to be 1
    
    for (t in 1:N){
      int betarow=1;
      for (i in 1:K){     // previous state
        for (j in 1:K){   // current state
          if (i==j){
            unA[t,i,j]=1; // exp(0)=1
          } else {
            unA[t,i,j]=exp(beta[betarow]*to_vector(u[t])); // exp(beta_0 + beta_1*sin + beta_2*cos)
            betarow=betarow+1;
          }
          
        }
        
      }
      for (i in 1:K)
        A[t][i]=unA[t][i]/sum(unA[t][i]);
        
    }

  // hyperparameters for mean and variance  
  for (k in 1:K){
    target+= gamma_lpdf(mu[k] | mu0, k0);
    target+= scaled_inv_chi_square_lpdf(sigma2y[k] | v0, sigma02);
    //target+= inv_gamma_lpdf(sigma2y[k] | v0, sigma02);
  }
  
  // priors for multinomial coefficients in tpm
  for (i in 1:K*(K-1)){
    for (j in 1:M)
        beta[i,j] ~ normal(0, 10);
        
  }
  
  // At each time point, there's a transition probability matrix. Each entry in a K by K matrix
  // is a linear combination of covariates (circadian components)
  {  
    real acc[K];
    real gamma[N, K];
   
    // forward algorithm
  
    // for (k in 1:K)
    // gamma[1, k] = normal_lpdf(y[1] | mu[k], sqrt(sigma2y[k]));
    // initial probabilities
    // for (k in 1:K) {
    //   pi1[k] ~ dirichlet(alpha);
    // }
    for (k in 1:K)
      gamma[1, k] = log(pi1[k]) + normal_lpdf(y[1] | mu[k], sqrt(sigma2y[k]));
    
    for (t in 2:N) {
      for (k in 1:K) {  // current state
        for (j in 1:K){ // previous state
        // accumulator
        // belief state + transition prob + local evidence at t
          acc[j] = gamma[t-1, j] + log(A[t][j,k]) + normal_lpdf(y[t] | mu[k], sqrt(sigma2y[k]));
        }
        gamma[t, k] = log_sum_exp(acc);
      }
    }
    target += log_sum_exp(gamma[N]); // log sum of posterior mixture likelihood
    
  }
  
}

