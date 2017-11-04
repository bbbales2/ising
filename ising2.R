library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
sigma = 0.01
S = 1000

source("ising_helpers2.R")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

ising_gibbs_derivs(x, 0.0, beta, S, 1)

# This is the output for the grid of test parameters
ising_sweep = function(x, beta, S) {
  list(mu = seq(-5.0, 5.0, length = 21),
       seed = 1:4) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      ising_gibbs_derivs(x, row$mu, beta, S, row$seed) %>%
        mutate(mu = row$mu, seed = row$seed, S = S) %>%
        mutate(!!!beta)
      }, mc.cores = 12) %>%
    bind_rows %>%
    mutate(seed = factor(seed))
}

samples = list()

for(i in 1:10000) {
  beta = rnorm(5, 0.0, 0.25)
  names(beta) = c("b0", "b1", "b2", "b3", "b4")
  samples[[i]] = ising_sweep(x, beta, S)
  print(i)
}

fit_iid = stan('models/iid.stan', data = list(N = nrow(data), y = data$phi, sigma = 0.1), iter = 10)
fit_mvn = stan('models/mvn.stan', data = list(N = nrow(data), y = data$phi, x = data$mu, sigma = 0.1, l = 0.5), iter = 10)

lpres = results %>% group_by(beta, gamma, seed) %>%
  summarize(nlp = -log_prob(fit_iid, mphi),
            dnlp_dbeta = -grad_log_prob(fit_iid, mphi) %*% dmphi_dbeta,
            dnlp_dgamma = -grad_log_prob(fit_iid, mphi) %*% dmphi_dgamma)

mlpres = lpres %>%
  group_by(beta, gamma) %>%
  summarize(nlp = mean(nlp))

lpres %>%
  mutate(angle = atan2(-dnlp_dgamma, -dnlp_dbeta),
         mag = sqrt(dnlp_dgamma^2 + dnlp_dbeta^2)) %>%
  #filter(log(nlp) < 7.5) %>%
  mutate(length = 0.025 * mag / max(mag),
         x = beta - cos(angle) * length / 2.0,
         y = gamma - sin(angle) * length / 2.0,
         xend = beta + cos(angle) * length / 2.0,
         yend = gamma + sin(angle) * length / 2.0) %>%
  ggplot(aes(x, y)) +
  geom_contour(data = mlpres, aes(x = beta, y = gamma, z = nlp, color = ..level..), bins = 100) +
  geom_segment(aes(xend = xend, yend = yend, alpha = 0.25),
               arrow = arrow(length = unit(0.01, "npc"))) +
  xlab("beta") + ylab("gamma") +
  geom_point(data = as.tibble(list(beta = betaRef, gamma = gammaRef)),
             aes(x = beta, y = gamma), color = "red", shape = "x", size = 5.0) +
  ggtitle(paste("Contours of negative log probability\nTruth at red x, beta = ", betaRef, ", gamma = ", gammaRef))

# Laplace approximation stuff
dataS %>% group_by(mu) %>% select(d2mphi_dbeta2) %>% summarize_all(funs(mean, sd))
loss = results %>% filter(beta == betaRef, gamma == gammaRef) %>%
  group_by(mu) %>% summarize_at(vars(everything()), mean) %>%
  select(-seed) %>%
  left_join(data, by = "mu") %>%
  mutate(nlp = (mphi - phi)^2,
         dnlp_dbeta = 2 * (mphi - phi) * dmphi_dbeta,
         dnlp_dgamma = 2 * (mphi - phi) * dmphi_dgamma,
         d2nlp_dbeta2 = 2 * dmphi_dbeta^2 + 2 * (mphi - phi) * d2mphi_dbeta2,
         d2nlp_dgamma2 = 2 * dmphi_dgamma^2 + 2 * (mphi - phi) * d2mphi_dgamma2,
         d2nlp_dbetadgamma = 2 * dmphi_dbeta * dmphi_dgamma + 2 * (mphi - phi) * d2mphi_dbetadgamma) %>%
  summarize(nlp = sum(nlp),
            dnlp_dbeta = sum(dnlp_dbeta),
            dnlp_dgamma = sum(dnlp_dgamma),
            d2nlp_dbeta2 = sum(d2nlp_dbeta2),
            d2nlp_dgamma2 = sum(d2nlp_dgamma2),
            d2nlp_dbetadgamma = sum(d2nlp_dbetadgamma))

#lpres = 
results %>% group_by(mu, beta, gamma) %>% select(mu, mphi) %>%
  summarize_at(vars(everything()), mean) %>%
  left_join(data, by = "mu") %>%
  mutate(nlp = (mphi - phi)^2) %>%
  group_by(beta, gamma) %>%
  summarize(nlp = sum(nlp)) %>%
  ungroup() %>%
  mutate(loss = 0.5 * loss$d2nlp_dgamma2 * (gamma - gammaRef)^2 +
           0.5 * loss$d2nlp_dbeta2 * (beta - betaRef)^2 +
           loss$d2nlp_dbetadgamma * (beta - betaRef) * (gamma - gammaRef)) %>%
  ggplot(aes(x = beta, y = gamma)) +
  geom_contour(aes(z = loss), bins = 100, color = "red") +
  geom_contour(aes(z = nlp), bins = 100) +
  xlab("beta") + ylab("gamma") +
  geom_point(data = as.tibble(list(beta = betaRef, gamma = gammaRef)),
             aes(x = beta, y = gamma), color = "red", shape = "x", size = 5.0) +
  ggtitle(paste("Contours of negative log probability\nTruth at red x, beta = ", betaRef, ", gamma = ", gammaRef))
