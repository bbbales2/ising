library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.25
beta = 0.15
kT = 1.0
sigma = 0.01
S = 100000

sourceCpp("ising.cpp")

dmu = 0.00001
((ising_kmc_deriv(x, mu + dmu, beta, kT, 1000, 0)$phi - ising_kmc_deriv(x, mu, beta, kT, 1000, 0)$phi) / dmu)
dbeta = 0.00001
((ising_kmc_deriv(x, mu, beta + dbeta, kT, 1000, 0)$phi - ising_kmc_deriv(x, mu, beta, kT, 1000, 0)$phi) / dbeta)
ising_kmc_deriv(x, mu, beta, kT, 1000, 0)

dmu = 0.00001
((ising_kmc_deriv(x, mu + dmu, beta, kT, 2, 0)$phi - ising_kmc_deriv(x, mu, beta, kT, 2, 0)$phi) / dmu)
dbeta = 0.00001
((ising_kmc_deriv(x, mu, beta + dbeta, kT, 2, 0)$phi - ising_kmc_deriv(x, mu, beta, kT, 2, 0)$phi) / dbeta)
ising_kmc_deriv(x, mu, beta, kT, 2, 0)

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
system.time(ising_kmc_deriv(x, mu, beta, kT, S, 0))

mus = seq(-4.0, 4.0, length = 20)
y = mclapply(mus, function(mu) { ising_kmc_deriv(x, mu, beta, kT, S, 1)$phi }, mc.cores = 24) %>% unlist() + rnorm(20) * sigma

as.tibble(list(mus = mus, y = y)) %>% ggplot(aes(mus, y)) +
  geom_point()

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

UgradU = function(q) {
  beta = q[1]
  
  out = mclapply(mus, function(mu) {
    tmp = ising_kmc_deriv(x, mu, beta, kT, S, 0)
    tmp["x"] = NULL
    tmp %>% as.tibble
  }, mc.cores = 24) %>% bind_rows
  
  logp = sum(-(out$phi - y)^2 / sigma^2)
  dlogpdbeta = sum(-2 * (out$phi - y) * out$dphidbeta / sigma^2)
  list(u = -logp,
       dudq = -dlogpdbeta)
}
beta = 0.25
dbeta = 0.0000001
(UgradU(beta + dbeta)$u - UgradU(beta)$u) / dbeta
UgradU(beta)

betas = seq(0.3, 0.4, length = 100)
mlogp = lapply(betas, function(beta) { UgradU(beta)$u }) %>% unlist()
as.tibble(list(beta = betas, mlogp = mlogp)) %>% ggplot(aes(beta, mlogp)) +
  geom_point()

source("radford.R")

bcs = rep(0, 50)
bcs[1] = 0.25
for(i in 2:length(bcs)) {
  bcs[i] = radford(UgradU, 1e-4, 50, bcs[i - 1])
  print(bcs[i])
  guess = mclapply(mus, function(mu) { ising_kmc_deriv(x, mu, bcs[i], kT, S, 0)$phi }, mc.cores = 24) %>% unlist()
  print(as.tibble(list(mus = mus, y = y, guess = guess)) %>% ggplot(aes(mus, y)) +
    geom_point() +
    geom_line(aes(mus, guess), colour = "red"))
}
bcs

