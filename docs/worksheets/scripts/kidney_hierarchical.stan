// Hierarchical model with common parameters for state rates across counties
// 

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of counties
  int<lower=0> Ns; // number of states
  vector[N] counts; // population size (COLUMN)
  int<lower = 0> y[N]; // 10 year instances of cancer (ROW)
  int<lower = 1> s[N]; // state index

}

// The parameters accepted by the model. 
parameters {
  real tau;
  vector[Ns] gamma;
}

// The transformed parameters accepted by the model.
transformed parameters {
  vector[Ns] state_theta = inv_logit(tau + gamma);
}

// The model to be estimated. We model the output
model {
  tau ~ normal(-10.0,0.5);
  gamma ~ normal(0.0,0.5);
  for (i in 1:N) {
    y[i] ~ poisson(10*counts[i]*state_theta[s[i]]);
  }
}

// Quantities generated conditional on the posterior
generated quantities {
  // posterior predictive values
  vector[N] y_rep;
  for (i in 1:N) {
    y_rep[i] = poisson_rng(10*counts[i]*state_theta[s[i]]);
  }
  
  // prior predictive values check
  //real<upper=1> theta_prior = exp(-normal_rng(10,1));
  //array[N] real y_rep_prior = poisson_rng(10*counts*theta_prior);
  
  // log likelihood for each y_i for LOO package
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y[i] | 10*counts[i]*state_theta[s[i]]);
  }
}

