library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.2
beta = 0.15
kT = 1.0

sourceCpp("ising.cpp")

S = 100000 - 1
x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
out = ising(x, mu, beta, kT, S)

getphi = function(mu, beta, tol = 0.1, chains = 4) {
  S = 10^4 - 1
  x = list()
  states = list()

  for(chain in 1:chains) {
    xT = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
    out = ising(xT, mu, beta, kT, S) #Warmup
    out = ising(out$x, mu, beta, kT, S)
    x[[chain]] = out$x
    states[[chain]] = as.tibble(out$states) %>% mutate(chain = chain)
  }
  
  for(f in 1:100) {
    pc = bind_rows(states) %>%
      group_by(chain) %>%
      summarize(m = mean(solo),
                v = var(solo),
                n = n())
    
    #wt = mean(pc %>% pull(v))
    #bt = var(pc %>% pull(m))
    
    #sig2h = (S - 1) * wt / S + bt
    #nut = 2 * (sig2h + bt / chains)^2 / wt
    #r = ((sig2h + bt / chains) / wt) * nut / (nut - 2)
    #print.data.frame(pc)
    #print(paste(sig2h, wt, bt, f, r))
    
    if(sd(pc %>% pull(m)) < tol) {
      print(f)
      break;
    }
    
    for(chain in 1:chains) {
      out = ising(x[[chain]], mu, beta, kT, S)
      x[[chain]] = out$x
      states[[chain]] = bind_rows(states[[chain]], as.tibble(out$states) %>% mutate(chain = chain))
    }
  }
  
  bind_rows(states)
}

getphi(0.2, 0.5, 1.0) %>%
  group_by(chain) %>%
  top_n(10000) %>%
  mutate(rn = row_number()) %>%
  gather(key, value, c(solo, pairs)) %>%
  ggplot(aes(rn, value, color = key)) +
  geom_line() +
  facet_grid(. ~ chain)

getEdE = function(mu, beta, r = 1.03) {
  getphi(mu, beta, r)
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

getEdE(mu, beta)

S = 1000
x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
ising(x, mu, beta, kT, S)
as.tibble(out$states) %>% mutate(rn = row_number()) %>%
  gather(key, value, c(solo, pairs)) %>% ggplot(aes(rn, value, color = key)) + geom_point()

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
