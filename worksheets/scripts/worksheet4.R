library(readr)
library(dplyr)
library(ggplot2)
library(usmap)
library(rstan)
library(bayesplot)
library(loo)

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

