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

N = 2 * 3 * 5 * 2
sigma = 0.01
S = 100
mus = seq(-5.0, 5.0, length = 21)

source("ising_helpers2.R")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

(ising_gibbs(x, 0.25, beta2, S, 2) -> out)$x %>% p

(y = ising_gibbs_derivs(x, 5.0, c(0.0, 0.0, 0.0, 0.0, 0.0), S, 1))

# This is the output for the grid of test parameters
ising_sweep = function(x, beta, S, seeds) {
  list(mu = mus,
       seed = seeds) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      ising_gibbs_derivs(x, row$mu, beta, S, row$seed) %>%
        mutate(mu = row$mu, seed = row$seed, S = S) %>%
        mutate(!!!beta)
    }, mc.cores = 24) %>%
    bind_rows %>%
    mutate(seed = factor(seed))
}

K = rbf_cov_vec(mus %>% as.matrix, mus %>% as.matrix, c(1.0))

beta2 = rnorm(5, 0.1, 0.25)
names(beta2) = c("b0", "b1", "b2", "b3", "b4")
data = list(X0 = ising_sweep(x, beta2, S * 100, 2) %>% pull(X0),
            Q = ising_gibbs_derivs(x, 0.0, beta2, S * 100, 2) %>%
              gather(name, Q, starts_with("Q")) %>% pull(Q))

fit_iid = stan('models/iid.stan',
               data = list(N = length(data$X0),
                           M = length(data$Q),
                           X0 = data$X0,
                           Q = data$Q,
                           useQ = 1,
                           sigma = 0.1),
               iter = 1)

opts = list()
for(o in 1:50) {
  b = rnorm(5, 0.1, 0.25)
  names(b) = c("b0", "b1", "b2", "b3", "b4")
  bs = list()
  us = list()
  
  for(i in 1:100) {
    u = UgradU(b)
    bs[[i]] = b
    us[[i]] = u$u
    
    cat(u$u, "|", b, "|", u$dudq, "\n")
    
    b = b - 0.01 * u$dudq / sqrt(sum(u$dudq^2))
  }
  
  opts[[o]] = do.call("cbind", bs) %>%
    t %>% as.tibble %>%
    mutate(which_opt = o, lp = -(us %>% unlist))
}

opt_df = opts %>% bind_rows

pmap(combn(names(b), 2) %>% t %>% as.tibble,
     ~ opt_df %>% select(..1, ..2, which_opt, lp) %>%
       rename(x = !!..1, y = !!..2) %>%
       mutate(which_x = !!..1,
              which_y = !!..2)) %>%
  bind_rows %>%
  ggplot(aes(x, y)) +
  geom_point(aes(group = which_opt, colour = lp), size = 0.01) +
  facet_grid(which_x ~ which_y)

bind_rows(list(x = mus, y = ising_sweep(x, sample(bs[50:100], 1) %>% unlist, S, sample(10000000, 1)) %>% pull(X0), which = "fit") %>% as.tibble,
          list(x = mus, y = data$X0, which = "data") %>% as.tibble) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(colour = which))
