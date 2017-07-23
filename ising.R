library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
require(Rcpp)

p = function(x) {
  as.tibble(melt(x)) %>%
    ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
    geom_tile() +
    xlab("x") +
    ylab("y") +
    coord_equal() +
    scale_y_reverse()
}

E1 = function(x) {
  tot = 0.0
  
  for(i in 1:nrow(x)) {
    for(j in 1:ncol(x)) {
      tot = tot + mu * x[i, j];
      
      if(i < nrow(x)) {
        tot = tot + beta * x[i, j] * x[i + 1, j]
      }
        
      if(j < ncol(x)) {
        tot = tot + beta * x[i, j] * x[i, j + 1]
      }
    }
  }
  
  -tot
}

solo = function(x) {
  sum(x)  
}

pairs = function(x) {
  tot = 0.0
  
  for(i in 1:nrow(x)) {
    for(j in 1:ncol(x)) {
      if(i < nrow(x)) {
        tot = tot + x[i, j] * x[i + 1, j]
      }
      
      if(j < ncol(x)) {
        tot = tot + x[i, j] * x[i, j + 1]
      }
    }
  }
  
  tot
}

E = function(x) {
  -mu * solo(x) - beta * pairs(x)
}

deltaE = function(x, i, j) {
  tot = -2.0 * mu * x[i, j]
  if(i < nrow(x)) {
    tot = tot + -2.0 * beta * x[i, j] * x[i + 1, j]
  }
  
  if(i > 1) {
    tot = tot + -2.0 * beta * x[i - 1, j] * x[i, j]
  }
  
  if(j < ncol(x)) {
    tot = tot + -2.0 * beta * x[i, j] * x[i, j + 1]
  }
  
  if(j > 1) {
    tot = tot + -2.0 * beta * x[i, j - 1] * x[i, j]
  }
    
  -tot
}

dsolo = function(x, i, j) {
  -2.0 * x[i, j]
}

dpairs = function(x, i, j) {
  tot = 0.0
  
  if(i < nrow(x)) {
    tot = tot + -2.0 * x[i, j] * x[i + 1, j]
  }
  
  if(i > 1) {
    tot = tot + -2.0 * x[i - 1, j] * x[i, j]
  }
  
  if(j < ncol(x)) {
    tot = tot + -2.0 * x[i, j] * x[i, j + 1]
  }
  
  if(j > 1) {
    tot = tot + -2.0 * x[i, j - 1] * x[i, j]
  }

  tot
}

deltaE = function(x, i, j) {
  -mu * dsolo(x, i, j) - beta * dpairs(x, i, j)
}

e = E(x)

for(i in 1:N) {
  for(j in 1:N) {
    y = x
    y[i, j] = y[i, j] * -1
    if(!all.equal((E(y) - E(x)), deltaE(x, i, j))) {
      print(paste(i, j))
      print(paste(E(y) - E(x), deltaE(x, i, j)))
    }
  }
}

N = 10
mu = 0.2
beta = 0.5
kT = 1.0

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

sourceCpp("ising.cpp")
e = E(x)
run = function(mu, beta, x, S) {
  out = array(0, dim = c(S, 2))
  colnames(out) = c("solo", "pairs")
  
  solo_ = solo(x)
  pairs_ = pairs(x)
  for(s in 1:S) {
    i = sample(N, 1)
    j = sample(N, 1)
    
    dsolo_ = dsolo(x, i, j)
    dpairs_ = dpairs(x, i, j)
    dE = mu * dsolo_ + beta * dpairs_
    r = runif(1)
    
    #print(paste(dE, exp(-dE / kT), r))
    if(dE < 0.0 | exp(-dE / kT) > r) {
      #print("hi")
      e = e + dE
      solo_ = solo_ + dsolo_
      pairs_ = pairs_ + dpairs_
      x[i, j] = x[i, j] * -1
    }
    
    out[s, 1] = solo_
    out[s, 2] = pairs_
  }
  
  list(x = x,
       states = out)
}
system.time(run(mu, beta, x, S))
system.time(ising(x, mu, beta, kT, S))
S = 1000000
out = ising(x, mu, beta, kT, S)
x = out$x
p(x)
as.tibble(out$states) %>%
  mutate(group = sample(S, S) %/% 101) %>%
  group_by(group) %>%
  summarize(mean_solo = mean(solo),
            mean_pairs = mean(pairs),
            sd_solo = sd(solo),
            sd_pairs = sd(pairs),
            conc = mean(solo / (2 * N * N) + 1.0 / 2.0),
            dmdmu = -cov(solo, solo) / kT,
            dmdbeta = -cov(solo, pairs) / kT,
            n = n()) %>%
  select(mean_solo, conc, dmdmu, dmdbeta) %>%
  summarize_all(funs(mean, sd))

as.tibble(out$states) %>%
  mutate(rn = row_number()) %>%
  sample_n(1000) %>%
  gather(type, number, c(solo, pairs)) %>%
  ggplot(aes(rn, number, color = type)) +
  geom_point()
p(x)

mus = seq(-4.0, 4.0, length = 100)

results = bind_rows(map(mus, function(mu) {
  out = run(mu, beta, x, S)
  x = out$x
  out = run(mu, beta, x, S)
  x = out$x
  as.tibble(out$states) %>%
    mutate(group = row_number() %/% 1001) %>%
    group_by(group) %>%
    summarize(mean_solo = mean(solo),
              mean_pairs = mean(pairs),
              sd_solo = sd(solo),
              sd_pairs = sd(pairs),
              conc = mean(solo / (2 * N * N) + 1.0 / 2.0),
              cov_solo = -cov(solo, solo) / kT,
              cov_solo_pairs = -cov(solo, pairs) / kT,
              n = n()) %>%
    select(conc) %>%
    summarize_all(funs(mean, sd))
}))

results %>% mutate(mu = mus) %>%
  ggplot(aes(mu, mean)) +
  geom_ribbon(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd)) +
  geom_line() + 
  coord_flip()


getEdE = function(mu, beta) {
  S = 10000000
  x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
  out = run(mu, beta, x, S)
  x = out$x
  out = run(mu, beta, x, S)
  x = out$x
  as.tibble(out$states) %>%
    mutate(group = sample(S, S) %/% 10001) %>%
    group_by(group) %>%
    summarize(mean_solo = mean(solo),
              mean_pairs = mean(pairs),
              sd_solo = sd(solo),
              sd_pairs = sd(pairs),
              conc = mean(solo / (2 * N * N) + 1.0 / 2.0),
              dmdmu = -cov(solo, solo) / kT,
              dmdbeta = -cov(solo, pairs) / kT,
              n = n()) %>%
    select(mean_solo, conc, dmdmu, dmdbeta) %>%
    summarize(m_phi = mean(mean_solo),
              m_dmdmu = mean(dmdmu),
              m_dmdbeta = mean(dmdbeta),
              sd_mphi = sd(mean_solo),
              sd_dmdmu = sd(dmdmu),
              sd_dmdbeta = sd(dmdbeta))
}

s1 = getEdE(0.2, 0.5)
s2 = getEdE(0.2, 0.5 + 0.01)
s1
s2
(s2 - s1) / 0.01
