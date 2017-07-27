library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.25
beta = 0.15
kT = 1.0

sourceCpp("ising.cpp")

getphi = function(mu, beta, kT, tol = 1.01, chains = 4) {
  S = 10^4 - 1
  x = list()
  states = list()
  
  for(chain in 1:chains) {
    xT = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
    out = ising_reduced(xT, mu, beta, kT, S) #Warmup
    out = ising_reduced(out$x, mu, beta, kT, S)
    x[[chain]] = out$x
    states[[chain]] = as.tibble(out$states) %>% mutate(chain = chain)
  }
  
  for(f in 1:100) {
    r = as.tibble(pc %>% summarize(mt = mean(m),
                  bt = sum((m - mt)^2) / (chains - 1),
                  wt = mean(v),
                  n = n[1]) %>%
      mutate(r = sqrt(((n - 1) * wt / n + bt) / wt)) %>%
      pull(r)
    
    if(r < tol) {
      break;
    }
    
    for(chain in 1:chains) {
      out = ising(x[[chain]], mu, beta, kT, S)
      x[[chain]] = out$x
      states[[chain]] = bind_rows(states[[chain]], as.tibble(out$states) %>% mutate(chain = chain))
    }
  }
  
  bind_rows(states) %>% summarize(u = mean(solo),
                                  grad = -cov(solo, pairs) / kT)
}

system.time(getphi(0.2, 0.5, 1.01))

getphi(0.2, 0.5, kT, 1.001)

#%>%
  group_by(chain) %>%
  top_n(10000) %>%
  mutate(rn = row_number()) %>%
  gather(key, value, c(solo, pairs)) %>%
  ggplot(aes(rn, value, color = key)) +
  geom_line() +
  facet_grid(. ~ chain)
