library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)
library(GGally)

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
beta = rnorm(5, 0.0, 0.1)

source("ising_helpers2.R")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

(ising_gibbs(x, 0.25, beta2, S, 2) -> out)$x %>% p

noise_df = map(1:100, ~ ising_gibbs_derivs(x, 0.0, beta2, S, .)) %>% bind_rows

noise_df %>% summarize_all(sd) %>% select(-starts_with("Q"))

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

names(beta) = c("b0", "b1", "b2", "b3", "b4")
sweep = ising_sweep(x, beta, S, 1:10)

K = rbf_cov_vec(mus %>% as.matrix, mus %>% as.matrix, c(1.0))

sweep[, 1:8] %>% filter(seed != 1) %>%
  gather(var, y, 1:6) %>%
  group_by(var, seed) %>%
  mutate(yh = K %*% fsolve(K + diag(0.25, length(mus)), y %>% as.matrix)) %>%
  ggplot(aes(mu, y)) +
  geom_point(aes(color = seed)) +
  geom_line(aes(y = yh, color = seed)) +
  facet_wrap(~ var, scales = "free_y")

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

b = rnorm(5, 0.1, 0.25)
names(b) = c("b0", "b1", "b2", "b3", "b4")
source('radford.R')

UgradU = function(b) {
  df = ising_sweep(x, b, S, sample(10000000, 1))
  
  grad = grad_log_prob(fit_iid, c(df$X0, b))
  gradv = grad[1:length(df$X0)]
  gradb = grad[(length(df$X0) + 1):length(grad)]
  
  list(u = -attr(grad, "log_prob"),
       dudq = -gradb -
         df %>% gather(var, y, starts_with("dX0")) %>%
           group_by(var) %>%
           mutate(yh = K %*% fsolve(K + diag(0.25, length(mus)), y %>% as.matrix)) %>%
           summarize(y = y %*% gradv,
                     yh = yh %*% gradv) %>% pull(yh))
}

bind_rows(list(x = mus, y = ising_sweep(x, sample(bs[50:100], 1) %>% unlist, S, sample(10000000, 1)) %>% pull(X0), which = "fit") %>% as.tibble,
          list(x = mus, y = data$X0, which = "data") %>% as.tibble) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(colour = which))

bs = list(b)
b[[2]] = radford(UgradU, 1e-3, 10, b)
b[[3]] = radford(UgradU, 1e-3, 10, b)
b[[4]] = radford(UgradU, 1e-3, 10, b)

for(i in 5:105) {
  b = radford(UgradU, 1e-3, 50, b)
  bs[[i]] = b
}

u = lapply(bs, UgradU)
list(s = 1:length(u),
     u = map(u, ~ .x$u) %>% unlist) %>% as.tibble %>%
  ggplot(aes(s, u)) +
  geom_point()

(posterior = do.call("cbind", bs) %>%
    t %>% as.tibble %>%
    mutate(lp = map(u, ~ .x$u) %>% unlist) %>%
    filter(row_number() > 50))

posterior %>%
  ggpairs(mapping = ggplot2::aes(fill = lp),
          lower = list(continuous = wrap("points", alpha = 0.3, size=0.1)))

posterior %>% filter(b0 > -1 & b0 < 1) %>% summary

for(i in 1:50) {
  u = UgradU(b)
  
  cat(u$u, "|", b, "|", u$dudq, "\n")
  
  b = b - 0.01 * u$dudq / sqrt(sum(u$dudq^2))
}
