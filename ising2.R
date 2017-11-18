library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)
library(GGally)

sourceCpp("covariance2.cpp")

p = function(x) {
  as.tibble(melt(x)) %>%
    ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
    geom_tile() +
    xlab("x") +
    ylab("y") +
    coord_equal() +
    scale_y_reverse()
}

N = factorial(5)
sigma = 0.01
S = 100
mus = seq(-5.0, 5.0, length = 21)
beta = rnorm(5, 0.0, 0.1)

source("ising_helpers2.R")

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

(ising_gibbs(x, 5.0, c(1.0, 1.0, 1.0, 1.0, 1.0), S, 1) -> out)$x %>% p

(y = ising_gibbs_derivs(x, 5.0, c(0.0, 0.0, 0.0, 0.0, 0.0), S, 1))

# This is the output for the grid of test parameters
ising_sweep = function(x, beta, S, seeds) {
  list(mu = mus,
       seed = seeds) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      ising_gibbs_derivs(x, row$mu, beta, S, row$seed) %>%
        mutate(mu = row$mu, seed = row$seed, S = S) %>%
        mutate(!!!beta)
    }, mc.cores = 24) %>%
    bind_rows %>%
    mutate(seed = factor(seed))
}

names(beta) = c("b0", "b1", "b2", "b3", "b4")
sweep = ising_sweep(x, beta, S, 1:10)

K = rbf_cov_vec(mus %>% as.matrix, mus %>% as.matrix, c(1.0))

sweep[, 1:8] %>% filter(seed != 1) %>%
  gather(var, y, 1:6) %>%
  group_by(var, seed) %>%
  mutate(yh = K %*% fsolve(K + diag(0.25, length(mus)), y %>% as.matrix)) %>%
  ggplot(aes(mu, y)) +
  geom_point(aes(color = seed)) +
  geom_line(aes(y = yh, color = seed)) +
  facet_wrap(~ var, scales = "free_y")

beta2 = rnorm(5, 0.1, 0.25)
names(beta2) = c("b0", "b1", "b2", "b3", "b4")
data = ising_sweep(x, beta2, S, 2) %>% select(-starts_with("dX0"))

fit_iid = stan('models/iid.stan', data = list(N = nrow(data), y = data$X0, sigma = 0.1), iter = 1)

b = rnorm(5, 0.1, 0.25)
names(b) = c("b0", "b1", "b2", "b3", "b4")
source('radford.R')

UgradU = function(b) {
  df = ising_sweep(x, b, S, sample(10000000, 1))
  
  grad = grad_log_prob(fit_iid, df$X0)
  
  list(u = -attr(grad, "log_prob"),
       dudq = -df %>% gather(var, y, starts_with("dX0")) %>%
       group_by(var) %>%
       mutate(yh = K %*% fsolve(K + diag(0.25, length(mus)), y %>% as.matrix)) %>%
       summarize(y = y %*% grad,
                 yh = yh %*% grad) %>% pull(yh))
}

(b = radford(UgradU, 1e-3, 50, b))

for(i in 1:10) {
  
  
  
  
  cat(attr(grad, "log_prob"), "|", b, "|", db, "\n")
  
  g = bind_rows(df %>% mutate(name = "current"),
                data %>% mutate(name = "data")) %>%
    ggplot(aes(mu, X0)) +
    geom_line(aes(color = name))
  print(g)
  Sys.sleep(0)
  
  b = b + 0.05 * db / sqrt(sum(db^2))
}

lp = log_prob(fit_iid, df$X0)
for(i in 1:5) {
  bn = b + rnorm(5, 0.0, 0.001)
  df = ising_sweep(x, bn, S, sample(10000000, 1))
  
  lpn = log_prob(fit_iid, df$X0)
  
  cat(lpn, "|", lp, "|", b, "|", bn, "\n")
  
  if(lpn > lp | runif(1) < exp(lpn - lp)) {
    b = bn;
    lp = lpn;
    cat("accept\n")
  } else {
    cat("reject\n")
  }
}

df %>% select()

(sweep %>% group_by(sample, seed) %>%
    summarize(lp = log_prob(fit_iid, X0),
              b0 = b0[1],
              b1 = b1[1],
              b2 = b2[1],
              b3 = b3[1],
              b4 = b4[1]) -> dfn) %>%
  ungroup() %>%
  select(-sample, -seed, -nlp) %>%
  summary()

ungroup() %>%
  select(-sample) %>%
  filter(nlp < 1.0) %>%
  ggpairs()

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
  select(-seed) %>%
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
