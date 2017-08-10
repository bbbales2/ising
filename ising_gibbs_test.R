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
S = 10000

sourceCpp("ising.cpp")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
ising(x, mu, beta, kT, S * N * N, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))
ising_gibbs(x, mu, beta, 0.0, kT, S * N * N, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))
