source("casm_helper.R")


source("ising_helpers.R")
require(Rcpp)
sourceCpp("likelihoods.cpp")

library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(stats)
library(parallel)

path = "/home/bbales2/casm/cubic_2d"
setECIs(path, c(0.0, 0.0, 1.0, 0.5, 1.0))
runMC(path)
corrs = getCorrs(path)
results = getResults(path)

corrs %>% ggplot(aes(mci, corr1)) +
  geom_point(aes(color = as.factor(t)))

results %>% ggplot(aes(param_chem_pota, corr1)) +
  geom_point()

out = list(beta = seq(-0.5, 1.5, length = 21),
     gamma = seq(-0.5, 1.5, length = 21)) %>%
  expand.grid %>%
  as.tibble %>%
  mutate(output = pmap(., function(beta, gamma) {
    setECIs(path, c(0.0, 0.0, beta, gamma))
    runMC(path)
    list(results = getResults(path), corrs = getCorrs(path))
  })) %>%
  mutate(results = map(output, ~ .x$results),
         corrs = map(output, ~ .x$corrs)) %>%
  select(-output)

out %>%
  select(-corrs) %>%
  unnest %>%
  group_by(beta, gamma) %>%
  select(param_chem_pota, corr1) %>%
  ggplot(aes(param_chem_pota, corr1)) +
  geom_point() +
  facet_grid(beta ~ gamma, labeller = label_both)

save(list = ls(all.names = TRUE), file = paste0(path, ".RData"), envir = .GlobalEnv)
