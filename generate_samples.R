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
for(i in 1:10000) {
  ecis = c(c(0.0, 0.0), rnorm(3))
  setECIs(path, ecis)
  runMC(path)
  corrs = getCorrs(path)
  results = getResults(path)
  saveRDS(list(ecis = ecis, corrs = corrs, results = results), paste0(path, "/samples/", i, ".rds"))
}
