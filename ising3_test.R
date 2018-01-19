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

N = 20
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
