library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)

N = 10
mu = 0.0
betas = seq(-0.5, 0.75, length = 21)
gammas = seq(-0.5, 0.75, length = 21)
kT = 1.0
sigma = 0.01
S = 100

sourceCpp("ising.cpp")

ising_gibbs(x, mu, 0.5, -0.25, kT, S, 0)$x
x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

# This is the output for the grid of test parameters
results = list(mu = seq(-5.0, 5.0, length = 21),
               beta = betas,
               gamma = gammas,
               seed = sample(1:1000, 4, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, beta, gamma, seed) {
    ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
      summarize(mphi = mean(solo) / N^2,
                mphi1 = mean(pairs) / N^2,
                dmphi_dbeta = -cov(solo, pairs) / kT / N^2,
                dmphi1_dbeta = -var(pairs) / kT / N^2,
                dmphi_dgamma = -cov(solo, corners) / kT / N^2,
                dmphi1_dgamma = -cov(pairs, corners) / kT / N^2) %>%
      mutate(mu = mu, beta = beta, gamma = gamma, seed = seed, S = S)
  }) %>% bind_rows %>%
  mutate(seed = factor(seed))

# This is the output for the reference parameters
betaRef = betas[which.min(abs(betas - 0.1))]
gammaRef = gammas[which.min(abs(gammas - 0.5))]
dataS = list(mu = results %>% pull(mu) %>% unique,
             seed = sample(1:1000, 4, replace = F)) %>%
  expand.grid %>%
  as.tibble %>%
  pmap(function(mu, seed) {
    ising_gibbs(x, mu, betaRef, gammaRef, kT, S * 100, seed)$states %>% as.tibble %>%
      summarize(phi = mean(solo) / N^2,
                phi1 = mean(pairs) / N^2,
                dmphi_dbeta = -cov(solo, pairs) / kT / N^2,
                dmphi1_dbeta = -var(pairs) / kT / N^2,
                dmphi_dgamma = -cov(solo, corners) / kT / N^2,
                dmphi1_dgamma = -cov(pairs, corners) / kT / N^2,
                d2mphi_dbeta2 = (mean(solo * pairs^2) +
                                 -mean(solo * pairs) * mean(pairs) +
                                 -cov(solo, pairs) * mean(pairs) +
                                 -var(pairs) * mean(solo)) / N^2,
                d2mphi_dgamma2 = (mean(solo * corners^2) +
                                  -mean(solo * corners) * mean(corners) +
                                  -cov(solo, corners) * mean(corners) +
                                  -var(corners) * mean(solo)) / N^2,
                d2mphi_dbetadgamma = (mean(solo * pairs * corners) +
                                      -mean(solo * pairs) * mean(corners) +
                                      -cov(solo, corners) * mean(pairs) +
                                      -cov(corners, pairs) * mean(solo)) / N^2) %>%
      mutate(mu = mu, seed = seed)
  }) %>% bind_rows

data = dataS %>% group_by(mu) %>% summarize_at(vars(everything()), mean) %>% select(-seed)

fit_iid = stan('models/iid.stan', data = list(N = nrow(data) + 1, y = c(data$phi, data$phi1[[10]]), sigma = 0.1), iter = 10)

lpres = results %>% group_by(beta, gamma, seed) %>%
  summarize(nlp = -log_prob(fit_iid, c(mphi, mphi1[[10]])),
            dnlp_dbeta = -grad_log_prob(fit_iid, c(mphi, mphi1[[10]])) %*% c(dmphi_dbeta, dmphi1_dbeta[[10]]),
            dnlp_dgamma = -grad_log_prob(fit_iid, c(mphi, mphi1[[10]])) %*% c(dmphi_dgamma, dmphi1_dgamma[[10]]))

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
  geom_segment(aes(xend = xend, yend = yend, alpha = 0.25),
               arrow = arrow(length = unit(0.01, "npc"))) +
  xlab("beta") + ylab("gamma") +
  geom_point(data = as.tibble(list(beta = betaRef, gamma = gammaRef)),
             aes(x = beta, y = gamma), color = "red", shape = "x", size = 5.0) +
  ggtitle(paste("Contours of negative log probability\nTruth at red x, beta = ", betaRef, ", gamma = ", gammaRef))

# Laplace approximation stuff
dataS %>% group_by(mu) %>% select(d2mphi_dbeta2) %>% summarize_all(funs(mean, sd))
loss = results %>% filter(beta == betaRef, gamma == gammaRef) %>%
  group_by(mu) %>% summarize_at(vars(everything()), mean) %>%
  select(-seed, -dmphi_dbeta, -dmphi_dgamma) %>%
  left_join(data, by = "mu") %>%
  mutate(nlp = (mphi - phi)^2,
         dnlp_dbeta = 2 * (mphi - phi) * dmphi_dbeta,
         dnlp_dgamma = 2 * (mphi - phi) * dmphi_dgamma,
         d2nlp_dbeta2 = 2 * dmphi_dbeta^2 + 2 * (mphi - phi) * d2mphi_dbeta2,
         d2nlp_dgamma2 = 2 * dmphi_dgamma^2 + 2 * (mphi - phi) * d2mphi_dgamma2,
         d2nlp_dbetadgamma = 2 * dmphi_dbeta * dmphi_dgamma + 2 * (mphi - phi) * d2mphi_dbetadgamma) %>%
  summarize(nlp = sum(nlp),
            dnlp_dbeta = sum(dnlp_dbeta),
            dnlp_dgamma = sum(dnlp_dgamma),
            d2nlp_dbeta2 = sum(d2nlp_dbeta2),
            d2nlp_dgamma2 = sum(d2nlp_dgamma2),
            d2nlp_dbetadgamma = sum(d2nlp_dbetadgamma))

#lpres = 
results %>% group_by(mu, beta, gamma) %>% select(mu, mphi) %>%
  summarize_at(vars(everything()), mean) %>%
  left_join(data, by = "mu") %>%
  mutate(nlp = (mphi - phi)^2) %>%
  group_by(beta, gamma) %>%
  summarize(nlp = sum(nlp)) %>%
  ungroup() %>%
  mutate(loss = 0.5 * loss$d2nlp_dgamma2 * (gamma - gammaRef)^2 +
           0.5 * loss$d2nlp_dbeta2 * (beta - betaRef)^2 +
           loss$d2nlp_dbetadgamma * (beta - betaRef) * (gamma - gammaRef)) %>%
  ggplot(aes(x = beta, y = gamma)) +
  geom_contour(aes(z = loss), bins = 100, color = "red") +
  geom_contour(aes(z = nlp), bins = 100) +
  xlab("beta") + ylab("gamma") +
  geom_point(data = as.tibble(list(beta = betaRef, gamma = gammaRef)),
             aes(x = beta, y = gamma), color = "red", shape = "x", size = 5.0) +
  ggtitle(paste("Contours of negative log probability\nTruth at red x, beta = ", betaRef, ", gamma = ", gammaRef))
