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
S = 100

sourceCpp("ising.cpp")

# This is the output for the grid of test parameters
results = list(mu = seq(-5.0, 5.0, length = 21),
               beta = seq(0.0, 1.0, length = 21),
               gamma = seq(0.0, 1.0, length = 21),
               seed = sample(1:1000, 4, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, beta, gamma, seed) {
    ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
      summarize(m_phi = mean(solo) / N^2,
                d_m_phi_d_beta = -cov(solo, pairs) / kT / N^2,
                d_m_phi_d_gamma = -cov(solo, corners) / kT / N^2) %>%
      mutate(mu = mu, beta = beta, gamma = gamma, seed = seed, S = S)
  }) %>% bind_rows %>%
  mutate(seed = factor(seed))

# This is the output for the reference parameters
beta = 0.1
gamma = 0.5
data = results %>% pull(mu) %>% unique %>%
  map(function(mu) {
    ising_gibbs(x, mu, beta, gamma, kT, S * 10, 0)$states %>% as.tibble %>%
      summarize(phi = mean(solo) / N^2) %>%
      mutate(mu = mu)
  }) %>% bind_rows

fit_iid = stan('models/iid.stan', data = list(N = nrow(data), y = data$phi, sigma = 0.1), iter = 10)
fit_mvn = stan('models/mvn.stan', data = list(N = nrow(data), y = data$phi, x = data$mu, sigma = 0.1, l = 0.5), iter = 10)
m_phi = results %>% filter(seed == 700 & beta == 0 & gamma == 0) %>% pull(m_phi)
grad_log_prob(fit_iid, m_phi)
grad_log_prob(fit_mvn, m_phi)

data %>% ggplot(aes(mu, phi)) +
  geom_point()

results %>% select(beta, gamma) %>% unique

lpres = results %>% group_by(beta, gamma, seed) %>%
  summarize(nlp = -log_prob(fit_mvn, m_phi),
            dnlp_dbeta = -grad_log_prob(fit_mvn, m_phi) %*% d_m_phi_d_beta,
            dnlp_dgamma = -grad_log_prob(fit_mvn, m_phi) %*% d_m_phi_d_gamma)

mlpres = lpres %>%
  group_by(beta, gamma) %>%
  summarize(nlp = mean(nlp))

lpres %>%
  mutate(angle = atan2(-dnlp_dgamma, -dnlp_dbeta),
         mag = sqrt(dnlp_dgamma^2 + dnlp_dbeta^2)) %>%
  #filter(log(nlp) < 7.5) %>%
  mutate(length = 0.025 * mag / max(mag),
         x = beta - cos(angle) * length / 2.0,
         y = gamma - sin(angle) * length / 2.0,
         xend = beta + cos(angle) * length / 2.0,
         yend = gamma + sin(angle) * length / 2.0) %>%
  ggplot(aes(x, y)) +
  geom_contour(data = mlpres, aes(x = beta, y = gamma, z = nlp, color = ..level..), bins = 100) +
  #scale_colour_gradientn(colours = rev(rainbow(4))) +
  #scale_colour_gradientn(colours = rev(terrain.colors(8))) +
  geom_segment(aes(xend = xend, yend = yend, alpha = 0.25),
               arrow = arrow(length = unit(0.01, "npc"))) +
  xlab("beta") + ylab("gamma") +
  geom_point(data = as.tibble(list(beta = beta, gamma = gamma)),
             aes(x = beta, y = gamma), color = "red", shape = "x", size = 5.0) +
  ggtitle(paste("Contours of negative log probability\nTruth at red x, beta = ", beta, ", gamma = ", gamma))

lpres %>%
  mutate(angle = atan2(-dnlp_dgamma, -dnlp_dbeta),
         mag = sqrt(dnlp_dgamma^2 + dnlp_dbeta^2)) %>%
  mutate(length = 0.01 * mag / max(mag),
         x = beta - cos(angle) * length / 2.0,
         y = gamma - sin(angle) * length / 2.0,
         xend = beta + cos(angle) * length / 2.0,
         yend = gamma + sin(angle) * length / 2.0) %>%
  arrange(length)

results

+ geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.3,"cm")))