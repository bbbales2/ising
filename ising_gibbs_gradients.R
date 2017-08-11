library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = -0.25
beta = 0.15
gamm = 0.25
kT = 1.0
sigma = 0.01
S = 10000

sourceCpp("ising.cpp")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)
ising(x, mu, beta, kT, S * N * N, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))
ising_gibbs(x, mu, beta, 0.0, kT, S * N * N, 0)$states %>% as.tibble %>%
  summarize(m_phi = mean(solo),
            sd_phi = sd(solo))

out = mclapply(sample(1:1000, 8, replace = F), function(s) {
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
  ggplot(aes(rn, value)) +
  geom_point(aes(color = which)) +
  facet_grid(~ seed)

results = list(mu = seq(-0.5, 0.5, length = 11),
               beta = seq(0.0, 0.25, length = 11),
               gamma = seq(0.0, 0.25, length = 11),
               seed = sample(1:1000, 3, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, beta, gamma, seed) {
    ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
      summarize(m_phi = mean(solo),
                d_m_phi_d_beta = -cov(solo, pairs) / kT,
                d_m_phi_d_gamma = -cov(solo, corners) / kT) %>%
      mutate(mu = mu, beta = beta, gamma = gamma, seed = seed)
    }) %>% bind_rows %>%
  mutate(seed = factor(seed))

beta = 0.10
gamma = 0.15
data = results %>% pull(mu) %>% unique %>%
  map(function(mu) {
    ising_gibbs(x, mu, beta, gamma, kT, S * 100, 0)$states %>% as.tibble %>%
      summarize(phi = mean(solo)) %>%
      mutate(mu = mu)
  }) %>% bind_rows

results %>% group_by(mu, beta, gamma) %>%
  summarize(sd_m_phi = sd(m_phi),
            m_phi = mean(m_phi),
            sd_d_m_phi_d_beta = sd(d_m_phi_d_beta),
            d_m_phi_d_beta = mean(d_m_phi_d_beta),
            sd_d_m_phi_d_gamma = sd(d_m_phi_d_gamma),
            d_m_phi_d_gamma = mean(d_m_phi_d_gamma))

lpres = results %>% left_join(data) %>%
  mutate(nlp = 2 * log(abs(m_phi - phi) + 1e-5),
         dnlp_dbeta = (2 / (abs(m_phi - phi) + 1e-5)) * d_m_phi_d_beta,
         dnlp_dgamma = (2 / (abs(m_phi - phi) + 1e-5)) * d_m_phi_d_gamma) %>%
  group_by(beta, gamma, seed) %>%
  summarize(nlp = sum(nlp),
            dnlp_dbeta = sum(dnlp_dbeta),
            dnlp_dgamma = sum(dnlp_dgamma))

mlpres = lpres %>%
  group_by(beta, gamma) %>%
  summarize(nlp = mean(nlp))


lpres %>%
  mutate(angle = atan2(dnlp_dgamma, dnlp_dbeta),
         length = 0.01,
         x = beta - cos(angle) * length / 2.0,
         y = gamma - sin(angle) * length / 2.0,
         xend = beta + cos(angle) * length / 2.0,
         yend = gamma + sin(angle) * length / 2.0,
         mag_grad = sqrt(dnlp_dgamma^2 + dnlp_dbeta^2)) %>%
  ggplot(aes(x, y)) +
  geom_tile(data = mlpres, aes(x = beta, y = gamma, fill = nlp)) +
  #scale_colour_gradientn(colours = rev(rainbow(4))) +
  scale_colour_gradientn(colours = rev(terrain.colors(8))) +
  geom_segment(aes(xend = xend, yend = yend, color = seed),
               arrow = arrow(length = unit(0.01, "npc"))) +
  xlab("beta") + ylab("gamma")

results

+ geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.3,"cm")))