library(MASS)

source("casm_helper.R")
source("ising_helpers.R")
require(Rcpp)

library(argparse)
library(reshape2)
library(tidyverse)
library(stats)
library(parallel)
library(gtools)
library(ggplot2)
library(rstan)

#there are n different pairs from which to build the triplet, so three opportunities
path = "/home/bbales2/casm/invent3"
N = 30
ecis = rep(0, length(getECIs(path)))
nonZero = c(2, 3, 4, 5, 6, 7, 14, 15)#, 16, 17, 18)
keep = c(2, 3, 4, 5, 6, 7)#, 14, 15, 16, 17, 18)
ts = c(1.0)
ecis[nonZero] = c(0.0, 0.440, 0.280, 0.100, -0.100, -0.133, 0.40, -0.1)#, 0.170, 0.080, -0.080)
mus = seq(-4, 4, 0.5)

paths = list()

for(i in 1:26) {
  paths[i] = paste0("/home/bbales2/ising/paper_outputs/individual4/ind", i, ".dat")
}

env = new.env()
load(paths[[1]], envir = env)

proc = map(paths, function(path) {
  env = new.env()
  load(path, envir = env)
  
  out = list(ecis = env$ecis[env$nonZero[-1]],
             ecis_est = env$b[env$nonZero[-1]])

  out
})

map(proc, ~ .$ecis - .$ecis_est) %>%
  do.call(rbind, .) %>%
  as.tibble %>%
  gather(key, value) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_grid(. ~ key) +
  ggtitle("Histogram of errors in ECI estimates (eci_true - eci_est)\nV6 and V7 are triplets, the rest are pairs\nUsing only pairs/single data")

map(proc, ~ sqrt(sum((.$ecis - .$ecis_est)^2))) %>%
  unlist %>%
  as.tibble %>%
  ggplot(aes(value)) +
  geom_histogram() +
  ggtitle("Histogram of sqrt of sum of squares difference\nbetween truth and estimate")

map(proc, ~ sqrt(sum(.$ecis^2))) %>%
  unlist %>%
  as.tibble %>%
  ggplot(aes(value)) +
  geom_histogram() +
  ggtitle("Histogram of lengths of eci vectors")

env$opts %>% map(~.[env$keep]) %>%
  do.call(rbind, .) %>%
  as.tibble %>%
  mutate(mu = mus) %>%
  gather(key, value, -mu) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = key))

map(data[[1]], ~ .$Eg[keep]) %>%
  do.call(rbind, .) %>%
  as.tibble %>%
  mutate(mu = mus) %>%
  gather(key, value, -mu) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = key))

list(error = map(proc, ~ sqrt(sum((.$ecis - .$ecis_est) ^ 2))) %>% unlist,
     triplets = map(proc, ~ sqrt(sum(.$ecis[6:7] ^ 2))) %>% unlist) %>%
  as.tibble %>%
  ggplot(aes(error, triplets)) +
  geom_point()
