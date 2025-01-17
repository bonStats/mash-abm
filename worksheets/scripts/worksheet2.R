library(readr)
library(ggplot2)
library(MCMCpack)
library(tidyr)

josh_data_dir <- function(fl){
  paste0("/Users/jbon/github/mash-abm/data/", fl)
}

golf <- read_csv(file = josh_data_dir("golf2.txt"))


golf %>% group_by(distance) %>% summarise(prop_success = mean(success)) %>%
  ggplot() + 
  geom_col(aes(x = distance, y = prop_success)) +
  scale_x_continuous("Distance to hole") +
  scale_y_continuous("Proportion successful") +
  theme_bw()

# prior specification
mean_beta0 <- 0
mean_beta1 <- -10
sd_beta0 <- 5
sd_beta1 <- 5

# beta 1 prior
ggplot(data.frame(x=c(-50,50)), aes(x=x)) + 
  stat_function(fun = dnorm, args=list(mean=mean_beta1, sd=sd_beta1)) + 
  theme_bw()

logit_samples <- MCMClogit(success ~ distance,
                                data = golf,
                                burnin = 1000,
                                mcmc = 10000,
                                thin = 10, # typically don't need to thin, can test what happens to acfplot() when thin = 1 or 10
                                verbose = 0,
                                b0 = c(mean_beta0, mean_beta1),
                                B0 = c(sd_beta0^(-2), sd_beta1^(-2))
)

# mean, variance, quantiles
summary(logit_samples)

# ESS
coda::effectiveSize(logit_samples)

# Autocorrelation plot
coda::acfplot(logit_samples)

# HPD Credible interval
HPDinterval(logit_samples)


## plot posterior logistic function

link <- function(x, b0, b1){
  exp(b0 + b1 * x) / ( 1 + exp(b0 + b1 * x) )
}

tidy_samples <- as_tibble(logit_samples) %>% 
  mutate(iter = 1:n()) %>%
  rename(beta0 = `(Intercept)`, beta1 = distance)

tidy_samples_theta <- expand_grid(
  distance = seq(0, max(golf$distance), length.out = 200), # x values to evaluate function
  iter = 1:nrow(tidy_samples) # iter index to join samples
  ) %>% 
  left_join(tidy_samples, by = "iter") %>% # join
  mutate(theta = link(distance, beta0, beta1))

# plot all posterior samples and use opacity (alpha) to visualise uncertainty
ggplot(tidy_samples_theta) + 
  geom_line(aes(x = distance, y = theta, group = iter), alpha = 0.01, colour = "blue") + 
  theme_bw()

# calculate pointwise 95% credible intervals of curve to visualise uncertainty
# include median
simdata %>% group_by(distance) %>%
  summarise(q025 = quantile(theta, 0.025), q975 = quantile(theta, 0.975), median = quantile(theta, 0.5)) %>%
  ggplot() +
  geom_line(aes(x = distance, y = q025), alpha = 1, colour = "blue") + 
  geom_line(aes(x = distance, y = q975), alpha = 1, colour = "blue") +
  geom_line(aes(x = distance, y = median), alpha = 1, colour = "red") +
  theme_bw()


# Multiple chains for R hat convergence diagnostic:
# Note: some R packages/functions will be able to run multiple chains at once for you.
# MCMClogit doesn't have this feature, so we will do it manually...

# The R hat convergence diagnostic uses
# "parallel chains ... with starting values that are overdispersed relative to the posterior distribution"
# If the prior is overdispersed relative to the posterior, use that!

logit_samples_list <- list()

for(r in 1:8){
  logit_samples_list[[r]] <- MCMClogit(success ~ distance,
                                  data = golf,
                                  burnin = 10000,
                                  mcmc = 15000,
                                  thin = 5,
                                  tune = 1.1,
                                  verbose = 0,
                                  beta.start = c(rnorm(1, mean = mean_beta0, sd = sd_beta0), rnorm(1, mean = mean_beta1, sd = sd_beta1)), #
                                  b0 = c(mean_beta0, mean_beta1),
                                  B0 = c(sd_beta0^(-2), sd_beta1^(-2))
  )
}


logit_samples2 <- as.mcmc.list(logit_samples_list)

# still too high
coda::gelman.diag(logit_samples2)

# try more samples

logit_samples_list <- list()

for(r in 1:8){
  logit_samples_list[[r]] <- MCMClogit(success ~ distance,
                                       data = golf,
                                       burnin = 20000,
                                       mcmc = 25000,
                                       thin = 5,
                                       tune = 1.1,
                                       verbose = 0,
                                       beta.start = c(rnorm(1, mean = mean_beta0, sd = sd_beta0), rnorm(1, mean = mean_beta1, sd = sd_beta1)), #
                                       b0 = c(mean_beta0, mean_beta1),
                                       B0 = c(sd_beta0^(-2), sd_beta1^(-2))
  )
}


logit_samples2 <- as.mcmc.list(logit_samples_list)

# it's ok!
coda::gelman.diag(logit_samples2)

gelman.plot(logit_samples2)

# we should redo our inference with these multiple chains because we are
# now more confident that they have converged

summary(logit_samples2)

# inference hasn't changed that much...


## Modelling from first principles


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


# bern log likelihood: y * log(prob) + (1-y) * log(1-prob)

log_likelihood1 <- function(sigma2,data,R) {
  
  p <- prob_success(sigma2, data$distance, R)
  sum(log(p) * data$success + log(1 - p) * (1 - data$success))
  
}

# more numerically stable
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


# To do: 
# Run multiple chains to calculate R hat to check convergence
# Check ESS, 
