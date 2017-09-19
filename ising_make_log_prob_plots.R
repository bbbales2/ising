library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = 0.0
betas = seq(-0.25, 0.75, length = 31)
gammas = seq(-0.25, 0.75, length = 31)
kT = 1.0
sigma = 0.01
S = 100

source("ising_helpers.R")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
# This is the output for the grid of test parameters
results = list(mu = seq(-5.0, 5.0, length = 21),
               beta = betas,
               gamma = gammas,
               seed = sample(1:1000, 4, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  (function(df) split(df, 1:nrow(df))) %>%
  mclapply(function(row) {
    ising_gibbs_derivs(x, row$mu, row$beta, row$gamma, kT, S, row$seed) %>%
      mutate(mu = row$mu, beta = row$beta, gamma = row$gamma, seed = row$seed, S = S)
  }, mc.cores = 24) %>%
  bind_rows %>%
  mutate(seed = factor(seed))

# This is the output for the reference parameters
mlpres = list(betaRef = betas[seq(5, 31, 7)],
              gammaRef = gammas[seq(5, 31, 7)]) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(betaRef, gammaRef) {
    dataS = (list(mu = results %>% pull(mu) %>% unique,
                  seed = sample(1:1000, 4, replace = F)) %>%
               expand.grid %>%
               as.tibble %>%
               (function(df) split(df, 1:nrow(df)))) %>%
      mclapply(function(row) {
        ising_gibbs_derivs(x, row$mu, betaRef, gammaRef, kT, S * 100, row$seed) %>% 
          mutate(mu = row$mu, seed = row$seed)
      }, mc.cores = 24) %>% bind_rows
    
    data = dataS %>%
      group_by(mu) %>%
      summarize_at(vars(everything()), mean) %>%
      rename(phi = mphi) %>%
      select(-seed) %>%
      select(mu, phi)
    
    fit_iid = stan('models/iid.stan', data = list(N = nrow(data), y = data$phi, sigma = 0.1), iter = 1)
    fit_mvn = stan('models/mvn.stan', data = list(N = nrow(data), y = data$phi, x = data$mu, sigma = 0.1, l = 0.5), iter = 1)
    
    lpres = results %>% group_by(beta, gamma, seed) %>%
      summarize(nlp = -log_prob(fit_iid, mphi),
                dnlp_dbeta = -grad_log_prob(fit_iid, mphi) %*% dmphi_dbeta,
                dnlp_dgamma = -grad_log_prob(fit_iid, mphi) %*% dmphi_dgamma)
    
    lpres %>%
      group_by(beta, gamma) %>%
      summarize(nlp = mean(nlp)) %>%
      mutate(betaRef = betaRef,
             gammaRef = gammaRef)
  }) %>% bind_rows %>% ungroup

mlpres %>%
  ggplot(aes(x = beta, y = gamma)) +
  geom_contour(aes(z = nlp, color = ..level..), bins = 100) +
  geom_point(data = mlpres %>% select(betaRef, gammaRef) %>% unique,
             aes(x = betaRef, y = gammaRef), color = "red", shape = "x", size = 5.0) +
  facet_grid(gammaRef ~ betaRef, labeller = label_both) +
  ggtitle(paste("Contours of negative log probability\nTruth at red x, beta = ", betaRef, ", gamma = ", gammaRef))
