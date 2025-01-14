library(readr)
library(dplyr)
library(ggplot2)
library(usmap)
library(rstan)
library(bayesplot)
library(loo)
library(brms)

josh_data_dir <- function(fl){
  paste0("/Users/jbon/github/mash-abm/data/", fl)
}

kidney <- read_csv(file = josh_data_dir("kidneycancerclean.csv"), skip = 4, col_types = cols(col_skip(), col_guess()))

# clean up
kidney <- kidney %>% 
  mutate(total_deaths = dc + dc.2) %>% # 1980-1989
  mutate(population = (pop+pop.2)/2) %>% # average in 1980-1989
  mutate(theta_hat = total_deaths/(10*population)) %>% # estimate of annual death rate
  mutate(low_rate = theta_hat <= quantile(theta_hat, 0.1)) %>%
  mutate(high_rate = theta_hat >= quantile(theta_hat, 0.9))

# plot example
plot_usmap("counties", data=kidney, values="high_rate") +
  scale_fill_discrete(h.start = 200, 
                      name = "High rate of kidney cancer deaths") 


# data summaries (using Frequentist estimate of annual rate of death)

kidney %>% 
  filter(theta_hat > 0 ) %>%
  mutate(logthetahat = log(theta_hat)) %>%
  summarise(logmean = mean(logthetahat), 
                     logmed = median(logthetahat), 
                     logsd = sd(logthetahat),
                     logq10 = quantile(logthetahat, 0.1), 
                     logq90 = quantile(logthetahat, 0.9))

scaling <- 10^5
kidney %>% group_by(state) %>% 
  summarise(mean = mean(theta_hat) * scaling, 
            med = median(theta_hat) * scaling, 
            q10 = quantile(theta_hat, 0.1) * scaling, 
            q90 = quantile(theta_hat, 0.9) * scaling,
            n = n()) 

## choose prior

# beta family

alpha <- 10^(-5)
beta <- 1

beta_mean <- alpha / (alpha + beta)
beta_sd <- sqrt((alpha * beta) / ( (alpha * beta + 1) * (alpha + beta)^2 ) )


## Stan fitting

kidney_stan_data = list(N = nrow(kidney),
                        counts = kidney$population,
                        y = kidney$total_deaths
)

fit_kidney = stan(file = "worksheets/scripts/kidney_global.stan", 
                  data = kidney_stan_data, iter = 4000, warmup = 3000,
                  )

print(fit_kidney, pars = c("theta","neg_log_theta"))

plot(fit_kidney, pars = c("neg_log_theta"))
plot(fit_kidney, pars = c("theta"))
pairs(fit_kidney, pars = c("neg_log_theta", "lp__"))

traceplot(fit_kidney, pars = c("neg_log_theta"))


# posterior predictive checks

y_rep <- as.matrix(fit_kidney, pars = "y_rep")

sample_id <- sample(1:kidney_stan_data$N, size = 8, replace = F)

ppc_hist(kidney_stan_data$y, y_rep[sample_id, ])
ppc_hist(log(1+kidney_stan_data$y), log(1+y_rep[sample_id, ]) )

# prior predictive checks

y_rep_prior <- as.matrix(fit_kidney, pars = "y_rep_prior")

sample_id <- sample(1:nrow(y_rep_prior), size = 8, replace = F)

ppc_hist(log(1+kidney_stan_data$y), log(1+y_rep_prior[sample_id, ]) )


# Bayes LOO 

loo(fit_kidney)


## Hierarchical model


kidney <- kidney %>%
  mutate(state_id = as.integer(factor(kidney$state)))
  

kidney_stan_data2 = list(N = nrow(kidney),
                        counts = kidney$population,
                        y = kidney$total_deaths,
                        s = kidney$state_id,
                        Ns = max(kidney$state_id)
)

fit_kidney = stan(file = "worksheets/scripts/kidney_hierarchical.stan", 
                  data = kidney_stan_data2, iter = 5000, warmup = 1000
)

 # needs more iterations! some gamma[i] are difficult to estimate


print(fit_kidney, pars = c("tau","gamma"))

pairs(fit_kidney, pars = c("tau", "gamma[31]"))


### Fit model with brms instead

## Global model

# Likelihood:
# y_i | theta ~ Pois(exp(Intercept + offset_i))
# offset_i = exp(offset) = 10 n_i (fixed)
# theta = Intercept (estimated)

# Prior check

prior_kidney_brms <- brm(total_deaths ~ 1 + offset(log(10*population)), 
                             family = poisson(link = "log"), 
                             sample_prior = "only", # NOT THE POSTERIOR
                             data = kidney)

prior_summary(prior_kidney_brms)

# Default prior:
# Intercept ~ student_t(df = 3, mean = -11.6026119141481, sd = 2.5)

prior_pred <- posterior_predict(prior_kidney_brms, ndraws = 10) # prior! but use the "posterior_predict()" function
ppc_hist(y = log(1+kidney$total_deaths), yrep = log(1+prior_pred))
ppc_dens_overlay(y = log(1+kidney$total_deaths), yrep = log(1+prior_pred))

log1_mean <- function(x) mean(log(1+x))
log1_sd <- function(x) sd(log(1+x))

pp_check(prior_kidney_brms, type = "stat", stat = "log1_mean")
pp_check(prior_kidney_brms, type = "stat", stat = "log1_sd")
pp_check(prior_kidney_brms, type = "stat_2d", stat = c("log1_mean","log1_sd"))


# Posterior fit

fit_kidney_brms <- brm(total_deaths ~ 1 + offset(log(10*population)), 
                       family = poisson(link = "log"), 
                       data = kidney)

summary(fit_kidney_brms)

# Posterior check

post_pred <- posterior_predict(fit_kidney_brms, ndraws = 8)
ppc_hist(y = log(1+kidney$total_deaths), yrep = log(1+post_pred))

ppc_dens_overlay(y = log(kidney$total_deaths), yrep = log(pred))
pp_check(fit_kidney_brms, type = "stat_2d")

pp_check(fit_kidney_brms, type = "stat", stat = "log1_mean")
pp_check(fit_kidney_brms, type = "stat", stat = "log1_sd")


## Hierarchical model

# Likelihood:
# y_i | theta ~ Pois(exp(Intercept + state_intercept[s[i]] + offset))
# offset = exp(offset) = 10 n_i (fixed)
# theta = Intercept (estimated)

# Prior check

prior_kidney_brms2 <- brm(total_deaths ~ (1|state) + offset(log(10*population)), 
                        family = poisson(link = "log"),
                        iter = 4000, warmup = 1000,
                        sample_prior = 'only',
                        data = kidney)

prior_summary(prior_kidney_brms2)

# Default prior
# Intercept ~ student_t(df = 3, mean = -11.6026119141481, sd = 2.5)
# state_intercept[s] = z[s] * sd[s]
# z[s] ~ N(0,1)
# sd[s] ~ student_t(3, 0, 2.5), Truncated > 0

prior_pred2 <- posterior_predict(prior_kidney_brms2, ndraws = 10) # prior! but use the "posterior_predict()" function
ppc_hist(y = log(1+kidney$total_deaths), yrep = log(1+prior_pred2))
ppc_dens_overlay(y = log(1+kidney$total_deaths), yrep = log(1+prior_pred2))

pp_check(prior_kidney_brms2, type = "stat", stat = "log1_mean")
pp_check(prior_kidney_brms2, type = "stat", stat = "log1_sd")
pp_check(prior_kidney_brms2, type = "stat_2d", stat = c("log1_mean","log1_sd"))

new_priors <- c(
  set_prior("normal(-11.6, 2)", class = "Intercept"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "state")
)

prior_kidney_brms2 <- brm(total_deaths ~ (1|state) + offset(log(10*population)), 
                          family = poisson(link = "log"),
                          iter = 4000, warmup = 1000,
                          sample_prior = 'only',
                          prior = new_priors,
                          data = kidney)

prior_pred2 <- posterior_predict(prior_kidney_brms2, ndraws = 19) # prior! but use the "posterior_predict()" function
ppc_hist(y = log(1+kidney$total_deaths), yrep = log(1+prior_pred2))
ppc_dens_overlay(y = log(1+kidney$total_deaths), yrep = log(1+prior_pred2))

pp_check(prior_kidney_brms2, type = "stat", stat = "log1_mean")
pp_check(prior_kidney_brms2, type = "stat", stat = "log1_sd")
pp_check(prior_kidney_brms2, type = "stat_2d", stat = c("log1_mean","log1_sd"))


# Posterior fit

fit_kidney_brms2 <- brm(total_deaths ~ (1|state) + offset(log(10*population)), 
                        family = poisson(link = "log"),
                        iter = 4000, warmup = 1000,
                        prior = new_priors,
                        data = kidney)

summary(fit_kidney_brms2) # model and main effect summary

ranef(fit_kidney_brms2) # random effect summary


# Posterior check

pp_check(fit_kidney_brms2, type = "hist", ndraws = 8)

post_pred2 <- posterior_predict(fit_kidney_brms2, ndraws = 19)
ppc_hist(y = log(1+kidney$total_deaths), yrep = log(1+post_pred2))

ppc_dens_overlay(y = log(1+kidney$total_deaths), yrep = log(1+post_pred2))

pp_check(fit_kidney_brms2, type = "stat_2d", stat = c("log1_mean","log1_sd"))

loo(fit_kidney_brms2, fit_kidney_brms)


# look at model fit by county

posterior_pred_mean <- colMeans(posterior_predict(fit_kidney_brms2))

kidney <- kidney %>% 
  mutate(pp_mean = posterior_pred_mean) %>%
  mutate(posterior_log_rate = log(pp_mean/(10*population)))

plot_usmap("counties", data=kidney, values="posterior_log_rate") +
  scale_fill_gradientn(colours = terrain.colors(3))
