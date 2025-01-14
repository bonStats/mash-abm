// Global rate model with common parameter for rates across counties
// 

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] counts; // population size (COLUMN)
  int<lower = 0> y[N]; // 10 year instances of cancer (ROW)
  
  // prior parameters for beta
  //real<lower=0> a;
  //real<lower=0> b;
}

// The parameters accepted by the model. 
parameters {
  real<lower=0> neg_log_theta;
}

// The transformed parameters accepted by the model.
transformed parameters {
  real<lower=0> theta = exp(-neg_log_theta);
}

// The model to be estimated. We model the output
model {
  neg_log_theta ~ normal(10,1); // choose prior
  y ~ poisson(10*counts*theta);
}

// Quantities generated conditional on the posterior
generated quantities {
  // posterior predictive values
  array[N] real y_rep = poisson_rng(10*counts*theta);
  
  // prior predictive values check
  real<upper=1> theta_prior = exp(-normal_rng(10,1));
  array[N] real y_rep_prior = poisson_rng(10*counts*theta_prior);
  
  // log likelihood for each y_i for LOO package
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y[i] | 10*counts[i]*theta);
  }
}

