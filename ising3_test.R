library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)
library(GGally)
require(gtools)

sourceCpp("covariance2.cpp")

p = function(x) {
  as.tibble(melt(x)) %>%
    ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
    geom_tile() +
    xlab("x") +
    ylab("y") +
    coord_equal() +
    scale_y_reverse()
}

N = 2 * 3 * 5
sigma = 0.01
S = 200
mus = seq(-5.0, 5.0, length = 21)

beta = c(0.0, 0.0, 0.0, 0.0, 0.0)
gamma = c(0.0)

source("ising_helpers3.R")
sourceCpp("ising3.cpp", showOutput = TRUE)

system.time(ising_gibbs(x, 0.0, beta, gamma, S, 0))

for(s in 1:100) {
  x = matrix(sample(c(-1, 1), size = N * N, replace = TRUE), nrow = N)
  
  i = sample(1:N, 1)
  j = sample(1:N, 1)
  
  n1 = triplets(x, 0)
  dn = dtriplets(x, i - 1, j - 1, 0)
  x[i, j] = -x[i, j]
  n2 = triplets(x, 0)
  
  stopifnot(n2 - n1 == -2 * dn)
  print(c(n2 - n1, -2 * dn))
}

ising_gibbs_derivs(x, 0.0, beta, gamma, S, 0)

ising_sweep = function(x, beta, gamma, S, seeds) {
  list(mu = mus,
       seed = seeds) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      res = ising_gibbs_derivs(x, row$mu, beta, gamma, S, row$seed)
      
      res$mu = row$mu
      res$seed = row$seed
      
      res
    }, mc.cores = 24)
}

{
  beta = rnorm(5, 0.1, 0.25)
  gamma = rnorm(1, 0.0, 0.25)
  data = map(ising_sweep(x, beta, gamma, S * 10, sample(1000000, 1)), ~ .$f)
  
  list(mus = mus,
       y = map(data, ~ .[[1]]) %>% unlist()) %>%
    as.tibble %>%
    ggplot(aes(mus, y)) +
    geom_line()
}
