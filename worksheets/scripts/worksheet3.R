golf

# |alpha - alpha_star| ~ Exp(lambda)
prob_success <- function(lambda,x,R){
  eps_tol <- atan(R/x)
  ( 1 - exp(-lambda*eps_tol)  ) / ( 1 - exp(-2*lambda*pi) )
}

ggplot(tibble(x=c(0.02,20)), aes(x=x)) + 
  stat_function(aes(colour = "Pr(Y=1)"), fun = prob_success, args = list(lambda = 100, R=0.0177)) + 
  scale_color_discrete("") + 
  theme_bw()

log_likelihood <- function(lambda,y,x,R) {
  
  sum(log(prob_success(lambda, x, R)) * y + log(1 - prob_success(lambda, x, R)) * (1 - y))
  
}


log_likelihood(27, y = golf$success, x = golf$distance, R = 0.177)


1/27


# z = (alpha - alpha_star)/(4pi) + 1/2 ~ Beta(gamma,gamma)
prob_success <- function(gamma,x,R){
  eps_tol <- atan(R/x)
  pbeta(0.5 + (eps_tol / (4*pi)), shape1 = gamma, shape2 = gamma) - pbeta(0.5 - (eps_tol / (4*pi)), shape1 = gamma, shape2 = gamma)
}

ggplot(tibble(x=c(0.02,20)), aes(x=x)) + 
  stat_function(aes(colour = "Pr(Y=1)"), fun = prob_success, args = list(gamma = 100000, R=0.0177)) + 
  scale_color_discrete("") + 
  theme_bw()

#  z = (alpha - alpha_star) ~ N(0,sigma2)
prob_success <- function(sigma,x,R){
  eps_tol <- atan(R/x)
  2*pnorm(eps_tol/sigma)  - 1
}

prob_success_logit <- function(x, b0 = 2.23, b1 = -0.25){
  exp(b0 + b1 * x) / ( 1 + exp(b0 + b1 * x) )
}

ggplot(tibble(x=c(0.02,20)), aes(x=x)) + 
  stat_function(aes(colour = "Pr(Y=1)"), fun = prob_success, args = list(sigma = 0.05, R = 0.177)) + 
  stat_function(aes(colour = "Pr(Y=1) logistic"), fun = prob_success_logit) +
  scale_color_discrete("") + 
  theme_bw()


library(MCMCpack)

# bern log likelihood: y * log(prob) + (1-y) * log(1-prob)

log_likelihood1 <- function(sigma2,data,R) {
  
  p <- prob_success(sigma2, data$distance, R)
  sum(log(p) * data$success + log(1 - p) * (1 - data$success))
  
}

log_likelihood2 <- function(sigma2,data,R) {
  
  p <- prob_success(sigma2, data$distance, R)
  prob_i <-ifelse(data$success==1, p,  1 - p)
  sum(log(prob_i))
  
}

log_likelihood1(0.1,golf,0.177)
log_likelihood2(0.1,golf,0.177)

log_likelihood1(0.01,golf,0.177)
log_likelihood2(0.01,golf,0.177)

log_prior <- function(sigma2) dgamma(sigma2, shape = 0.1, rate = 1, log = T)

log_post <- function(sigma2,data,R){

  log_likelihood2(sigma2,data,R) + log_prior(sigma2)

}


samples <- MCMCmetrop1R(fun = log_post, 
                        data = golf, 
                        R = 0.177, 
                        theta.init = 0.1, 
                        burnin = 5000, mcmc = 10000)


# to do: run multiple chains to calculate R hat to check convergence
# check ESS, s

