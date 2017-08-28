library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.0
beta = 0.15
gamma = 0.25
kT = 1.0
sigma = 0.01
S = 1000

sourceCpp("ising.cpp")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

ising(x, mu, beta, kT, S * N * N, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))
ising_gibbs(x, mu, beta, 0.0, kT, S, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))

beta = 0.15
gamma = 0.15
ising_gibbs(x, mu, beta, gamma, kT, S, 0)$x %>%
  melt %>%
  as.tibble %>%
  ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
  geom_tile() +
  xlab("x") +
  ylab("y") +
  coord_equal() +
  scale_y_reverse()

out = mclapply(sample(1:1000, 4, replace = F), function(s) {
  ising_gibbs(x, mu, beta, gamm, kT, S, s)$states %>% as.tibble %>%
    mutate(seed = s, rn = row_number())
}, mc.cores = 2) %>% bind_rows

out %>% group_by(seed) %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo),
            d_m_phi_d_beta = -cov(solo, pairs) / kT,
            d_m_phi_d_gamma = -cov(solo, corners) / kT) %>%
  summarize(sd_m_phi = sd(m_phi),
            m_phi = mean(m_phi),
            sd_d_m_phi_d_beta = sd(d_m_phi_d_beta),
            d_m_phi_d_beta = mean(d_m_phi_d_beta),
            sd_d_m_phi_d_gamma = sd(d_m_phi_d_gamma),
            d_m_phi_d_gamma = mean(d_m_phi_d_gamma))

out %>% gather(which, value, c(solo, pairs, corners)) %>%
  ggplot(aes(rn, value / N^2)) +
  geom_point(aes(color = which)) +
  facet_grid(~ seed)

results = list(mu = seq(-5.0, 5.0, length = 21),
               beta = c(-0.25, 0.0, 0.25, 0.5),
               gamma = c(-0.25, 0.0, 0.25, 0.5),
               seed = sample(1:1000, 4, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, beta, gamma, seed) {
    ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
      summarize(m_phi = mean(solo) / N^2) %>%
      mutate(mu = mu, beta = beta, gamma = gamma, seed = seed, S = S)
  }) %>% bind_rows %>%
  mutate(seed = factor(seed))

results %>% ggplot(aes(mu, m_phi)) +
  geom_point(aes(group = seed, color = seed), size = 1.0) +
  geom_line(aes(group = seed, color = seed)) +
  ylab("mean of phi") +
  facet_grid(beta ~ gamma, labeller = "label_both") +
  coord_flip()
