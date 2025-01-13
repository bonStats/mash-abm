//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] n; // population size
  int<lower = 0> y[N]; // 10 year instances of cancer
}

// The parameters accepted by the model. 
parameters {
  real<lower=0> theta;
}

// The transformed parameters accepted by the model.
transformed parameters {

}

// The model to be estimated. We model the output
model {
  theta ~ gamma(20, 430000);
  y ~ poisson(10*n*theta);
}

