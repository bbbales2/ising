library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.2
beta = 0.0
kT = 1.0

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

sourceCpp("ising.cpp")

getEdE = function(mu, beta) {
  S = 100000 - 1
  out = ising(x, mu, beta, kT, S)
  out = ising(out$x, mu, beta, kT, S)
  x = out$x
  as.tibble(out$states) %>%
    summarize(mean_solo = mean(solo),
              mean_pairs = mean(pairs),
              sd_solo = sd(solo),
              sd_pairs = sd(pairs),
              conc = mean(solo / (2 * N * N) + 1.0 / 2.0),
              dmdmu = -cov(solo, solo) / kT,
              dmdbeta = -cov(solo, pairs) / kT,
              n = n()) %>%
    select(mean_solo, conc, dmdmu, dmdbeta) %>%
    rename(mphi = mean_solo) %>%
    mutate(mu = mu)
}

s1 = getEdE(0.2, 0.5)
s2 = getEdE(0.2, 0.5 + 0.01)
s1
s2
(s2 - s1) / 0.01

mus = seq(-4.0, 4.0, length = 100)

results = bind_rows(mclapply(rep(mus, 10), function(mu) { getEdE(mu, beta) }, mc.cores = 24))
#summarize(m_phi = mean(mean_solo),
#          m_dmdmu = mean(dmdmu),
#          m_dmdbeta = mean(dmdbeta),
#          sd_mphi = sd(mean_solo),
#          sd_mdmdmu = sd(dmdmu),
#          sd_mdmdbeta = sd(dmdbeta))

results

results %>%
  gather(which, value, c(mphi, conc, dmdmu, dmdbeta)) %>%
  group_by(which, mu) %>%
  summarize(m = mean(value),
            l = min(value),
            u = max(value)) %>%
  ggplot(aes(mu, m)) +
  geom_ribbon(aes(ymin = l, ymax = u)) +
  geom_line(col = "red2") +
  facet_wrap(~ which, scales = "free_y")
