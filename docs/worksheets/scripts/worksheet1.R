library(readr)
library(dplyr)
library(ggplot2)

josh_data_dir <- function(fl){
  paste0("/Users/jbon/github/mash-abm/data/", fl)
}

deputes <- read_csv(file = josh_data_dir("deputes2019.csv"))
# change location to match folder of 'deputes2019.csv' on your computer
# e.g. "~/Downloads/deputes2019.csv" or "C:/Users/[User Name]/Downloads/deputes2019.csv"

ggplot(deputes) + 
  geom_histogram(aes(x = questions_orales)) +
  scale_x_continuous("Number of questions") +
  scale_y_continuous("Number of politicians") +
  #facet_wrap(~sexe) +
  theme_bw()



# Model 1: Gamma prior

alpha0 <- 3
beta0 <- 1

ybar <- mean(deputes$questions_orales)
n <- nrow(deputes)

# posterior
m1_alpha <- n * ybar + alpha0
m1_beta <- n + beta0


ggplot(data=tibble(lambda =c(0.001,5)), aes(x=lambda))+
  stat_function(fun=dgamma, args=list(shape=alpha0, rate=beta0), colour = "red") + # prior
  stat_function(fun=dgamma, args=list(shape=m1_alpha, rate=m1_beta), colour = "green") + # posterior
  theme_bw()

# mean
m1_alpha/m1_beta

# median
qgamma(0.5, shape = m1_alpha, rate = m1_beta)

# Credible intervals
c(qgamma(0.025, shape = m1_alpha, rate = m1_beta),
qgamma(1-0.025, shape = m1_alpha, rate = m1_beta))



# Model 2: Jeffery's prior

ybar1 <- mean(filter(deputes,sexe=="H")$questions_orales) # H
ybar2 <- mean(filter(deputes,sexe=="F")$questions_orales) # F


# posterior
m2_alpha1 <- n * ybar1 + alpha0
m2_alpha2 <- n * ybar2 + alpha0
m2_beta1 <- n + beta0
m2_beta2 <- n + beta0

ggplot(data=tibble(lambda =c(0.001,5)), aes(x=lambda))+
  stat_function(fun=dgamma, args=list(shape=alpha0, rate=beta0), colour = "red") + # prior
  stat_function(fun=dgamma, args=list(shape=m2_alpha1, rate=m2_beta1), colour = "green") + # posterior men
  stat_function(fun=dgamma, args=list(shape=m2_alpha2, rate=m2_beta2), colour = "orange") + # posterior women
  theme_bw()

# Credible intervals (Men)
c(qgamma(0.025, shape = m2_alpha1, rate = m2_beta1),
  qgamma(1-0.025, shape = m2_alpha1, rate = m2_beta1))

# Credible intervals (Women)
c(qgamma(0.025, shape = m2_alpha2, rate = m2_beta2),
  qgamma(1-0.025, shape = m2_alpha2, rate = m2_beta2))


# Model 1: Jeffery's prior

ybar <- mean(deputes$questions_orales)
n <- nrow(deputes)

# posterior
m1_alpha <- n * ybar + 0.5
m1_beta <- n


# Model 2: Jeffery's prior

ybar1 <- mean(filter(deputes,sexe=="H")$questions_orales) # H
ybar2 <- mean(filter(deputes,sexe=="F")$questions_orales) # F


# posterior
m2_alpha1 <- n * ybar1 + 0.5
m2_alpha2 <- n * ybar2 + 0.5
m2_beta1 <- n
m2_beta2 <- n


# Monte carlo approximation

# r_lambda = lambda_1 / lambda_2

N <- 100

# posterior
samples_lambda_1 <- rgamma(N, m2_alpha1, rate = m2_beta1)
samples_lambda_2 <- rgamma(N, m2_alpha2, rate = m2_beta2)

r_lambda <- samples_lambda_1 / samples_lambda_2
mean(r_lambda)


hist(r_lambda)

mean(r_lambda)
sqrt(var(r_lambda))

# prior 
samples_lambda_1 <- rgamma(N, alpha0, rate = beta0)
samples_lambda_2 <- rgamma(N, alpha0, rate = beta0)

r_lambda <- samples_lambda_1 / samples_lambda_2

hist(r_lambda)

mean(r_lambda)
sqrt(var(r_lambda))


