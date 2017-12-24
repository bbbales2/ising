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

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

ising_sweep = function(x, beta, S, seeds) {
  list(mu = mus,
       seed = seeds) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      res = ising_gibbs_derivs(x, row$mu, beta, S, row$seed)
      
      res$mu = row$mu
      res$seed = row$seed
      
      res
    }, mc.cores = 24)
}

data = map(ising_sweep(x, beta, S * 10, 2), ~ list(f = .$f))

UgradU = function(q) {
  u = 0.0
  dudq = rep(0, length(q))
  r = sample(10000000, 1)
  
  calc = ising_sweep(x, q, S, r)

  for(i in 1:length(calc)) {
    if(useFuncLoss) {
      u = u + 0.5 * (calc[[i]]$f[1] - data[[i]]$f[1])^2
      dudq = dudq + (calc[[i]]$f[1] - data[[i]]$f[1]) * calc[[i]]$jac[1,]
    }
    
    if(useCorrLoss) {
      u = u + 0.5 * sum((calc[[i]]$f[-1] - data[[i]]$f[-1])^2)
      dudq = dudq + (calc[[i]]$f[-1] - data[[i]]$f[-1]) %*% calc[[i]]$jac[-1,]
    }
  }
  
  list(u = 100 * u,
       dudq = 100 * dudq)
}

fn = function(b) { UgradU(b)$u }
gn = function(b) { UgradU(b)$dudq }

opts = list()
for(o in 1:100) {
  b = rnorm(5, 0.1, 0.25)
  names(b) = c("b0", "b1", "b2", "b3", "b4")
  out = optim(b, fn, gn, method = "L-BFGS-B", control = list(maxit = 200, trace = 1, REPORT = 1))
  opts[[o]] = out$par %>% as.list %>% as.tibble %>%
    mutate(lp = -out$value, which_opt = o)
  
  cat("finished optimization", o, "\n")
}

opt_samples = map(1:length(opts),
                  ~ ising_sweep(x, opts[[.]] %>% gather(name, b, 1:5) %>% pull(b) %>%
                                                 setNames(names(b)), S, sample(10000000, 1)))

save.image(args$output)