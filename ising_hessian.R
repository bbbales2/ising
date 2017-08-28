library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.0
beta = 0.15
gamma = 0.0
kT = 1.0
sigma = 0.01
S = 1000

sourceCpp("ising.cpp")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

ising_gibbs(x, mu, beta, 0.0, kT, S, 1)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))

ising_gibbs(x, mu, 0.5, 0.5, kT, S, 1)$x %>% melt %>%
  ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
  geom_tile() +
  xlab("x") +
  ylab("y") +
  coord_equal() +
  scale_y_reverse()

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
    mutate(seed = s, beta = beta)
}, mc.cores = 2) %>% bind_rows

full %>% ggplot(aes(beta, mphi)) +
  geom_line(aes(color = seed))

out %>% group_by(seed) %>%
  summarize(mphi = mean(solo),
            dmphi_dbeta = -cov(solo, pairs),
            dmphi_dgamma = -cov(solo, corners),
            d2mphi_dbeta2 = mean(solo * pairs^2) +
              -mean(solo * pairs) * mean(pairs) +
              -cov(solo, pairs) * mean(pairs) +
              -var(pairs) * mean(solo)) %>%
  select(-seed) %>%
  summarize_all(funs(mean, sd))

out %>% gather(which, value, c(solo, pairs, corners)) %>%
  ggplot(aes(rn, value)) +
  geom_point(aes(color = which)) +
  facet_grid(~ seed)

results = list(mu = seq(-0.5, 0.5, length = 11),
            beta = seq(0.10, 0.15, length = 21),
            gamma = seq(0.10, 0.15, length = 21),
            seed = sample(1:1000, 3, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, beta, gamma, seed) {
    ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
      summarize(mphi = mean(solo)) %>%
      mutate(mu = mu, beta = beta, gamma = gamma, seed = seed)
  }) %>% bind_rows %>%
  mutate(seed = factor(seed))

betaRef = results %>% pull(beta) %>% unique %>% median#0.4125#results %>% pull(beta) %>% unique %>% median
gammaRef = results %>% pull(beta) %>% unique %>% median#0.2250#results %>% pull(gamma) %>% unique %>% min
data = results %>% pull(mu) %>% unique %>%
  map(function(mu) {
    ising_gibbs(x, mu, betaRef, gammaRef, kT, S * 1000, 0)$states %>% as.tibble %>%
      summarize(phi = mean(solo),
                dmphi_dbeta = -cov(solo, pairs) / kT,
                dmphi_dgamma = -cov(solo, corners) / kT,
                d2mphi_dbeta2 = mean(solo * pairs^2) +
                  -mean(solo * pairs) * mean(pairs) +
                  -cov(solo, pairs) * mean(pairs) +
                  -var(pairs) * mean(solo),
                d2mphi_dbetadgamma = mean(solo * pairs * corners) +
                  -mean(solo * pairs) * mean(solo) +
                  -cov(solo, corners) * mean(pairs) +
                  -cov(solo, pairs) * mean(corners),
                d2mphi_dgamma2 = mean(solo * corners^2) +
                  -mean(solo * corners) * mean(corners) +
                  -cov(solo, corners) * mean(corners) +
                  -var(corners) * mean(solo)) %>%
      mutate(mu = mu)
  }) %>% bind_rows

lpres = results %>% left_join(data) %>%
  mutate(nlp = (mphi - phi)^2) %>%
  group_by(beta, gamma, seed) %>%
  summarize(nlp = sum(nlp))

lpres = results %>% left_join(data) %>%
  mutate(nlp = (mphi - phi)^2,
         dnlp_dbeta = 2 * (mphi - phi) * dmphi_dbeta,
         dnlp_dgamma = 2 * (mphi - phi) * dmphi_dgamma,
         d2nlp_dbeta2 = 2 * dmphi_dbeta^2 + 2 * (mphi - phi) * d2mphi_dbeta2,
         d2nlp_dgamma2 = 2 * dmphi_dgamma^2 + 2 * (mphi - phi) * d2mphi_dgamma2,
         d2nlp_dbetadgamma = 2 * dmphi_dbeta * dmphi_dgamma + 2 * (mphi - phi) * d2mphi_dbetadgamma) %>%
  group_by(beta, gamma, seed) %>%
  summarize(nlp = sum(nlp),
            dnlp_dbeta = sum(dnlp_dbeta),
            dnlp_dbeta = sum(dnlp_dgamma),
            d2nlp_dbeta2 = sum(d2nlp_dbeta2),
            d2nlp_dgamma2 = sum(d2nlp_dgamma2),
            d2nlp_dbetadgamma = sum(d2nlp_dbetadgamma))

mlpres = lpres %>%
  group_by(beta, gamma) %>%
  summarize(nlp = mean(nlp),
            dnlp_dbeta = mean(dnlp_dbeta),
            d2nlp_dbeta2 = mean(d2nlp_dbeta2),
            d2nlp_dgamma2 = mean(d2nlp_dgamma2),
            d2nlp_dbetadgamma = mean(d2nlp_dbetadgamma))

nlpres = lpres %>%
  filter((beta - betaRef)^2 < 1e-6 & (gamma - gammaRef)^2 < 1e-6) %>%
  summarise_all(mean)

nlpres

d2db2 = nlpres %>% pull(d2nlp_dbeta2) %>% mean
d2dg2 = nlpres %>% pull(d2nlp_dgamma2) %>% mean
d2dbg = nlpres %>% pull(d2nlp_dbetadgamma) %>% mean

mlpres %>% filter(beta > 0.0 & gamma > 0.0) %>%
  mutate(loss = 0.5 * d2dg2 * (gamma - gammaRef)^2 +
           0.5 * d2db2 * (beta - betaRef)^2 +
           d2dbg * (beta - betaRef) * (gamma - gammaRef)) %>%
  gather(type, value, c(nlp, loss)) %>%
  ggplot() +
  geom_contour(aes(x = beta, y = gamma, z = value), bins = 20) +
  xlab("beta") + ylab("gamma") +
  coord_equal() +
  facet_grid(. ~ type)
  
lpres %>% filter(gamma == gammaRef) %>%
  ggplot(aes(beta, nlp)) +
  geom_line(aes(color = seed)) +
  stat_function(fun = function(x) { (lpres %>% filter(beta == betaRef) %>% pull(d2nlp_dbeta2) %>% mean) * (x - betaRef)^2 })
