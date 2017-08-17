library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.0
beta = 0.15
gamm = 0.25
kT = 1.0
sigma = 0.01
S = 1000

sourceCpp("ising.cpp")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

ising_gibbs(x, mu, beta, 0.0, kT, S, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))

beta = 0.0
gamma = 0.5
ising_gibbs(x, mu, beta, gamma, kT, S, 0)$x %>%
  melt %>%
  as.tibble %>%
  ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
  geom_tile() +
  xlab("x") +
  ylab("y") +
  coord_equal() +
  scale_y_reverse()

out = mclapply(sample(1:1000, 8, replace = F), function(s) {
  ising_gibbs(x, mu, beta, gamma, kT, S, s)$states %>% as.tibble %>%
    mutate(seed = s, rn = row_number())
}, mc.cores = 2) %>% bind_rows

out %>% group_by(seed) %>%
  summarize(mphi = mean(solo),
            dmphi_dbeta = -cov(solo, pairs) / kT,
            dmphi_dgamma = -cov(solo, corners) / kT,
            d2mphi_dbeta2 = mean(solo * solo * pairs) -
              mean(solo * solo) * mean(pairs) -
              2 * mean(solo) * mean(solo * pairs) +
              2 * mean(solo)^2 * mean(pairs)) %>%
  select(-seed) %>%
  summarize_all(funs(mean, sd))

out %>% gather(which, value, c(solo, pairs, corners)) %>%
  ggplot(aes(rn, value)) +
  geom_point(aes(color = which)) +
  facet_grid(~ seed)

gamma = 0.15
full = list(mu = seq(-0.5, 0.5, length = 11),
            beta = seq(-0.25, 0.25, length = 11),
            seed = sample(1:1000, 3, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, beta, seed) {
    ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
      summarize(mphi = mean(solo)) %>%
      mutate(mu = mu, beta = beta, gamma = gamma, seed = seed)
  }) %>% bind_rows %>%
  mutate(seed = factor(seed))

beta = 0.0
data = results %>% pull(mu) %>% unique %>%
  map(function(mu) {
    ising_gibbs(x, mu, beta, gamma, kT, S, 0)$states %>% as.tibble %>%
      summarize(phi = mean(solo),
                dmphi_dbeta = -cov(solo, pairs) / kT,
                dmphi_dgamma = -cov(solo, corners) / kT,
                d2mphi_dbeta2 = mean(solo * solo * pairs) -
                  mean(solo * solo) * mean(pairs) -
                  2 * mean(solo) * mean(solo * pairs) +
                  2 * mean(solo)^2 * mean(pairs)) %>%
      mutate(mu = mu)
  }) %>% bind_rows

lpres = full %>% left_join(data) %>%
  mutate(nlp = (m_phi - phi)^2) %>%
  group_by(beta, seed) %>%
  summarize(nlp = sum(nlp))

full %>% filter(beta == 0.0) %>% left_join(data) %>%
  mutate(nlp = (mphi - phi)^2,
         dnlp_dbeta = 2 * (mphi - phi) * dmphi_dbeta,
         dnlp_dbeta = 2 * (mphi - phi) * dmphi_dbeta)

mlpres = lpres %>%
  group_by(beta) %>%
  summarize(nlp = mean(nlp))

lpres %>% ggplot(aes(beta, nlp)) +
  geom_line(aes(color = seed)) +
  stat_function(fun = function(x) {x^2})
