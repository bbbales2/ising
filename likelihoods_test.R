library(tidyverse)
library(rstan)
library(Rcpp)
library(Rcpp11)

sourceCpp("likelihoods.cpp")

squared_loss(c(1.0, 1.0, 3.0), seq(-1.0, 1.0, length = 3))
