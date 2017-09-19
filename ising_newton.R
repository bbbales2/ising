library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(stats)
library(parallel)
require(Rcpp)

N = 10
mu = 0.0
kT = 1.0
sigma = 0.01
S = 100

source("ising_helpers.R")
sourceCpp("likelihoods.cpp")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
# This is the output for the grid of test parameters
betaRef = 0.45
gammaRef = 0.15
data = list(mu = seq(-5.0, 5.0, length = 21)) %>%
  expand.grid %>%
  as.tibble %>%
  (function(df) split(df, 1:nrow(df))) %>%
  mclapply(function(row) {
    ising_gibbs_derivs(x, row$mu, betaRef, gammaRef, kT, S, 0) %>%
      mutate(mu = row$mu, beta = betaRef, gamma = gammaRef, S = S) %>%
      rename(phi = mphi) %>%
      select(mu, phi)
  }, mc.cores = 24) %>%
  bind_rows

func = function(beta, gamma) {
  sol = list(mu = seq(-5.0, 5.0, length = 21)) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      ising_gibbs_derivs(x, row$mu, beta, gamma, kT, S, 0) %>%
        mutate(mu = row$mu, beta = beta, gamma = gamma, S = S)
    }, mc.cores = 1) %>%
    bind_rows
  
  losses = squared_loss(data$phi, sol$mphi)
  
  sol %>% summarise(nlp = losses$nlp,
    dnlp_dbeta = losses$jacobian %*% dmphi_dbeta,
    dnlp_dgamma = losses$jacobian %*% dmphi_dgamma,
    d2nlp_dbeta2 = t(dmphi_dbeta) %*% losses$hessian %*% dmphi_dbeta + losses$jacobian %*% d2mphi_dbeta2,
    d2nlp_dgamma2 = t(dmphi_dgamma) %*% losses$hessian %*% dmphi_dgamma + losses$jacobian %*% d2mphi_dgamma2,
    d2nlp_dbetadgamma = t(dmphi_dbeta) %*% losses$hessian %*% dmphi_dgamma + losses$jacobian %*% d2mphi_dbetadgamma)
}

I = 10
mins = mclapply(1:12, function(r) {
  p = runif(2, -0.25, 0.75)
  betas = rep(0, I + 1)
  gammas = rep(0, I + 1)
  for(i in 1:I) {
    betas[i] = p[1]
    gammas[i] = p[2]
    UgradU = func(p[1], p[2])
    jac = c(UgradU$dnlp_dbeta, UgradU$dnlp_dgamma)
    hessian = matrix(c(UgradU$d2nlp_dbeta2, UgradU$d2nlp_dbetadgamma, UgradU$d2nlp_dbetadgamma, UgradU$d2nlp_dgamma2), nrow = 2)
    print(cat("nlogp: ", UgradU %>% pull(nlp), "p: ", p))
    p = p + tryCatch(solve(hessian, jac), error = function(e) NULL)
  }
  betas[I + 1] = p[1]
  gammas[I + 1] = p[2]
  
  list(beta = betas, gamma = gammas) %>% as.tibble %>% mutate(r = r, t = row_number())
}, mc.cores = 24) %>% bind_rows %>% mutate(r = as.factor(r)) %>%
  filter(beta > -0.25 & gamma > -0.25 & beta < 0.75 & gamma < 0.75)

logp = list(beta = seq(-0.25, 0.75, length = 21),
               gamma = seq(-0.25, 0.75, length = 21)) %>%
  expand.grid %>%
  as.tibble %>%
  (function(df) split(df, 1:nrow(df))) %>%
  mclapply(function(row) {
    func(row$beta, row$gamma) %>%
      mutate(beta = row$beta, gamma = row$gamma)
  }, mc.cores = 24) %>%
  bind_rows

logp %>% ggplot(aes(beta, gamma)) +
  geom_contour(aes(z = nlp), bins = 100, color = "black") +
  geom_line(data = mins, aes(group = r)) +
  geom_point(data = mins, aes(group = r, color = t), size = 1.0) +
  geom_point(data = as.tibble(list(beta = betaRef, gamma = gammaRef)),
           aes(x = beta, y = gamma), color = "red", shape = "x", size = 5.0) +
  ggtitle("Newton iteration to find minima for 12 minimizings\n(t is the iteration number)")
