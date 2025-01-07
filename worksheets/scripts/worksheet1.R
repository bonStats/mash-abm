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

y <- deputes$questions_orales
ybar <- mean(y)
n <- length(y)

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

# Model evidence 1

val_lambda <- 3.5
post_logdens1 <- function(lambda) dgamma(lambda, shape = m1_alpha, rate = m1_beta, log = T)
prior_logdens1 <- function(lambda) dgamma(lambda, shape = alpha0, rate = beta0, log = T)
loglikelihood1 <- function(lambda, y) sum(dpois(y, lambda = lambda, log = T)) # Don't use for MC, too small value
loglikelihoodC1 <- function(lambda, y) sum(y*log(lambda) - lambda)


log_model_evidence1 <-function(y, include_factorial = T){ # includes factorial term from likelihood by default
  
  ll <- ifelse(include_factorial, 
               loglikelihood1(val_lambda, y) , 
               loglikelihoodC1(val_lambda, y)
               )
  prior_logdens1(val_lambda) + ll - post_logdens1(val_lambda)
  
}

log_model_evidence1(deputes$questions_orales)
log_model_evidence1(deputes$questions_orales, include_factorial = F)

# Model 2: Gamma prior

y1 <- filter(deputes,sexe=="H")$questions_orales
y2 <- filter(deputes,sexe=="F")$questions_orales
ybar1 <- mean(y1) # H
ybar2 <- mean(y2) # F
n1 <- length(y1)
n2 <- length(y2)

# posterior
m2_alpha1 <- n1 * ybar1 + alpha0
m2_alpha2 <- n2 * ybar2 + alpha0
m2_beta1 <- n1 + beta0
m2_beta2 <- n2 + beta0

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


# Model evidence 1

val_lambda1 <- 3.5
val_lambda2 <- 3
post_logdens2 <- function(lambda1,lambda2)dgamma(lambda1, shape = m2_alpha1, rate = m2_beta1, log = T) + dgamma(lambda2, shape = m2_alpha2, rate = m2_beta2, log = T)
prior_logdens2 <- function(lambda1,lambda2) dgamma(lambda1, shape = alpha0, rate = beta0, log = T) + dgamma(lambda2, shape = alpha0, rate = beta0, log = T)
loglikelihood2 <- function(lambda1,lambda2, y1,y2) sum(dpois(y1,lambda1,log=T)) + sum(dpois(y2,lambda2,log=T))
loglikelihoodC2 <- function(lambda1,lambda2, y1,y2) sum(y1*log(lambda1) - lambda1) + sum(y2*log(lambda2) - lambda2)

log_model_evidence2 <- function(y1,y2, include_factorial = T){
  
  ll <- ifelse(include_factorial, 
               loglikelihood2(val_lambda1, val_lambda2, y1, y2), 
               loglikelihoodC2(val_lambda1, val_lambda2, y1, y2)
               )
  prior_logdens2(val_lambda1, val_lambda2) + ll - post_logdens2(val_lambda1, val_lambda2)
  
}

# BF: M2 vs M1

log_model_evidence2(y1, y2, include_factorial = T) -
  log_model_evidence1(y, include_factorial = T)

# log evidence excluding product of factorials from likelihood
log_model_evidence2(y1, y2, include_factorial = F) -
  log_model_evidence1(y, include_factorial = F)






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


