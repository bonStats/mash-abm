library(readr)
library(dplyr)
library(ggplot2)
library(usmap)
library(rstan)


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



## Stan

kidney_stan_data = list(N = nrow(kidney),
                        n = kidney$population,
                        y = kidney$total_deaths
)

fit_kidney = stan(file = "worksheets/scripts/kidney_iid.stan", data = kidney_stan_data)

fit_kidney

plot(fit_kidney)
pairs(fit_kidney, pars = c("theta", "lp__"))
traceplot(fit_kidney)
