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

source("ising_helpers2.R")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

(ising_gibbs(x, 0.25, beta2, S, 2) -> out)$x %>% p

ising_gibbs(x, 0.75, beta2, S, 2)$states %>% as.tibble %>%
  mutate(rn = row_number()) %>%
  gather(var, y, -rn) %>%
  ggplot(aes(rn, y)) +
  geom_point() +
  facet_grid(. ~ var)

noise_df = map(1:100, ~ ising_gibbs_derivs(x, 0.0, beta2, S, .)) %>% bind_rows

noise_df %>% summarize_all(sd) %>% select(starts_with("Q"))

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
data = list(X0 = ising_sweep(x, beta2, S * 10, 2) %>% pull(X0),
            Q = ising_gibbs_derivs(x, 0.0, beta2, S * 10, 2) %>%
              gather(name, Q, starts_with("Q")) %>% pull(Q))

fit_iid = stan('models/iid.stan',
               data = list(N = length(data$X0),
                           M = length(data$Q),
                           X0 = data$X0,
                           Q = data$Q,
                           useQ = 1,
                           sigma = 0.1),
               iter = 1)

UgradU = function(b) {
  r = sample(10000000, 1)
  df = ising_sweep(x, b, S, r)
  Q = ising_gibbs_derivs(x, 0.0, b, S, r) %>%
    gather(name, Q, starts_with("Q")) %>% pull(Q)
  
  grad = grad_log_prob(fit_iid, c(df$X0, Q))
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

opts = list()
for(o in 1:1) {
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
  
  cat("finished optimization", o, "\n")
}

opt_df = opts %>% bind_rows %>%
  group_by(which_opt) %>%
  top_n(1, lp) %>%
  mutate(type = "opt") %>%
  ungroup %>%
  bind_rows(beta2 %>% as.list %>% as.tibble %>% mutate %>%
    mutate(which_opt = 0, type = "truth", lp = 0, rn = 1)) %>%
  bind_rows(rnorm(5 * 2000, 0.1, 0.25) %>% matrix(ncol = 5) %>% as.tibble %>%
    setNames(names(b)) %>% mutate(which_opt = 0, type = "prior", lp = 0, rn = 1))

pairs = function(vec, n) {
  rbind(combn(names(b), 2) %>% t,
        rep(names(b), 2) %>% matrix(nrow = 5))
}

opt_plot = pmap(combn(names(b), 2) %>% t %>% as.tibble,
     ~ opt_df %>% select(..1, ..2, which_opt, lp, rn, type) %>%
       rename(x = !!..1, y = !!..2) %>%
       mutate(which_x = !!..1,
              which_y = !!..2)) %>%
  bind_rows

opt_dens = map(names(b) %>% as.list, ~ opt_df %>%
      rename(x = !!.x) %>%
      mutate(which = !!.x)) %>%
  bind_rows

opt_dens %>%
  ggplot(aes(x)) +
  geom_density(aes(colour = type, fill = type), alpha = 0.15) +
  facet_grid(. ~ which)

opt_plot %>% filter(type == "opt", lp > -10) %>%
  ggplot(aes(x, y)) +
  geom_density2d(data = opt_plot %>% filter(type == "prior"), bins = 10) +
  geom_point(aes(group = which_opt, colour = lp), size = 0.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  geom_point(data = opt_plot %>% filter(type == "truth"), colour = "black", size = 2.0, shape = 4, stroke = 2) +
  facet_grid(which_x ~ which_y)

opt_samples = map(1:length(opts), ~ list(x = mus,
                                         y = ising_sweep(x, opts[[.]] %>% top_n(1) %>%
                                                            gather(name, b, 1:5) %>% pull(b) %>%
                                                            setNames(names(b)), S, sample(10000000, 1)) %>% pull(X0),
                                         which = "opt", opt = ., lp = opts[[.]] %>% pull(lp) %>% max) %>% as.tibble) %>%
  bind_rows()

opt_samples %>% filter(lp > -10) %>%
  group_by(x) %>% summarize(mu = mean(y),
                            sd = sd(y)) %>%
  ggplot(aes(x, mu)) +
  geom_ribbon(aes(ymin = mu - 2 * sd, ymax = mu + 2 * sd), alpha = 0.5) +
  geom_line() +
  geom_point(data = list(x = mus, mu = data$X0, which = "data") %>% as.tibble, color = "red") +
  xlab("mu") + ylab("avg. composition") +
  ggtitle("truth are red dots,\ndistribution of data generated from fits in black/grey\nmean + 95% interval estimates")

x %>% p

ising_gibbs(x, 0.0, beta2, S, sample(10000000, 1))$x %>% p

ising_gibbs(x, 0.0, sample(opts, 1) %>% first %>% top_n(1) %>%
                    gather(name, b, 1:5) %>% pull(b) %>%
                    setNames(names(b)), S, sample(10000000, 1))$x %>% p
