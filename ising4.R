#!/usr/bin/env Rscript
library(argparse)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)
library(GGally)
require(gtools)

parser = ArgumentParser()

parser$add_argument("-f", action="store_true", default=FALSE, help="Use function loss")
parser$add_argument("-c", action="store_true", default=FALSE, help="Use correlation loss")
parser$add_argument("--config", default="", help = "File to be sourced that has all config stuff in it")
parser$add_argument("--output", default="", help = "Place to store the output workspace")

args = parser$parse_args()

if(length(args$config) == 0) {
  print("No configuration file provided")
  quit()
}

if(length(args$output) == 0) {
  print("No output filename provided")
  quit()
}

source(args$config)

useFuncLoss = args$f
useCorrLoss = args$c

if(!useFuncLoss & !useCorrLoss) {
  print("No loss function enabled")
  quit()
}

sourceCpp("covariance2.cpp")
source("ising_helpers2.R")

N = 2 * 3 * 5
sigma = 0.01
S = 200
mus = seq(-5.0, 5.0, length = 21)

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

# This is the output for the grid of test parameters
ising_sweep = function(x, beta, S, seeds) {
  list(mu = mus,
       seed = seeds) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      res = ising_gibbs_derivs(x, row$mu, beta, S, row$seed)
      df = res$jac[1,]
      names(df) = map(names(df), ~ paste0("dX0d", .))
      
      bind_cols(res$f %>% t %>% as.tibble, df %>% t %>% as.tibble) %>%
        mutate(mu = row$mu, seed = row$seed, S = S) %>%
        mutate(!!!beta)
    }, mc.cores = 24) %>%
    bind_rows %>%
    mutate(seed = factor(seed))
}

K = rbf_cov_vec(mus %>% as.matrix, mus %>% as.matrix, c(1.0))

data = list(X0 = ising_sweep(x, beta, S * 10, 2) %>% pull(X0),
            Q = ising_gibbs_derivs(x, 0.0, beta, S * 10, 2)$f[-1, ])

UgradU = function(b) {
  u = 0.0
  r = sample(10000000, 1)
  gradv = rep(0, length(b))
  gradb = rep(0, length(b))
  dudq = rep(0, length(b))
  
  if(useFuncLoss) {
    df = ising_sweep(x, b, S, r)
    u = u + 0.5 * sum((df$X0 - data$X0)^2)
    gradv = (df$X0 - data$X0)
    dudq = dudq + df %>% gather(var, y, starts_with("dX0")) %>%
      group_by(var) %>%
      mutate(yh = K %*% fsolve(K + diag(0.25, length(mus)), y %>% as.matrix)) %>%
      summarize(y = y %*% gradv,
                yh = yh %*% gradv) %>% pull(yh)
  }
  
  if(useCorrLoss) {
    dQ = ising_gibbs_derivs(x, 0.0, b, S, r)
    QQ = dQ$f[2:length(dQ$f), ]
    u = u + 0.5 * sum((QQ - data$Q)^2)
    dudq = dudq + (QQ - data$Q) %*% dQ$jac[2:length(dQ$f), ]
  }
  
  list(u = 100 * u,
       dudq = 100 * dudq)
}

fn = function(b) { UgradU(b)$u }
gn = function(b) { UgradU(b)$dudq }

opts = list()
for(o in 1:250) {
  b = rnorm(5, 0.1, 0.25)
  names(b) = c("b0", "b1", "b2", "b3", "b4")
  out = optim(b, fn, gn, method = "L-BFGS-B", control = list(maxit = 200, trace = 1, REPORT = 1))
  opts[[o]] = out$par %>% as.list %>% as.tibble %>%
    mutate(lp = -out$value, which_opt = o)
  
  cat("finished optimization", o, "\n")
}

opt_samples = map(1:length(opts), ~ bind_cols(list(x = mus,
                                                   which = "opt",
                                                   opt = .,
                                                   lp = opts[[.]] %>% pull(lp) %>% max) %>% as.tibble,
                                              ising_sweep(x, opts[[.]] %>% top_n(1) %>%
                                                            gather(name, b, 1:5) %>% pull(b) %>%
                                                            setNames(names(b)), S, sample(10000000, 1)))) %>%
  bind_rows()

save.image(args$output)
