library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

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

res = ising_gibbs_derivs(x, 0.0, beta2, S, 1)

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

tmp = ising_sweep(x, beta, S, 1:1)

data = do.call("cbind", map(ising_sweep(x, beta, S * 10, 2), ~ .$f))

map(ising_sweep(x, beta, S * 10, 2), ~ .$jac)

useCorrLoss = TRUE
useFuncLoss = TRUE
UgradU2 = function(b) {
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
  
  list(u = u,
       dudq = dudq)
}