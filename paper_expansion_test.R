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

path = "/home/bbales2/casm/invent1"
N = 30
ecis = rep(0, length(getECIs(path)))
nonZero = c(3, 4, 5, 6, 7, 14)
keep = c(3, 4, 5, 6, 7, 14)
ts = c(2.0)
ecis[nonZero[1]] = rnorm(1, 0.30, 0.1)
ecis[nonZero[2:length(nonZero)]] = rnorm(length(nonZero) - 1, 0.1, 0.25)
mus = seq(-4, 4, 0.5)

makeECIs = function() {
  ecis[keep] = rnorm(length(keep), 0.1, 0.25)
  ecis[-keep] = 0
  ecis
}

smu = function(corrs, Tfrac) {
  NN = N * N
  res = corrs %>% select(starts_with("corr")) %>% as.matrix
  
  jac = -cov(res, res) * NN / Tfrac
  rownames(jac) = colnames(res)
  colnames(jac) = colnames(res)
  f = colMeans(res)
  names(f) = colnames(res)
  
  out = list(g = t(res), Eg = f, Egrad = jac)
  
  out
}

runSimulation = function(path, g, mu = NULL) {
  setSupercell(path, N)
  setECIs(path, g)
  
  if(!is.null(mu)) {
    setChemicalPotentials(path, mu, mu, 0)
  } else {
    setChemicalPotentials(path, -4, 4, 0.5)
  }
  runMC(path)
  corrs = getCorrs(path)
  
  Tfrac = getTemperatureFraction(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500) %>%
    do(out = smu(., Tfrac)) %>% pull(out)
}

# Stolen from https://gist.github.com/doobwa/941125
log_sum_exp = function(x) {
  offset = max(x)
  log(sum(exp(x - offset))) + offset
}

#lambda = 0.0005
GgradG2 = function(path, g, data, mu = NULL) {
  a = runSimulation(path, g, mu)
  
  Tfrac = getTemperatureFraction(path)
  
  out = list()
  if(!is.null(mu)) {
    idxs = which(mus == mu)
  } else {
    idxs = 1:length(a)
  }
  
  for(i in idxs) {
    grad = rep(0, length(keep))
    grad2 = rep(0, length(keep))
    
    if(!is.null(mu)) {
      ai = 1
    } else {
      ai = i
    }
    
    out[[ai]] = list()
    
    out[[ai]]$lp = (-data[[i]]$Eg[keep] %*% g[keep] / Tfrac - log_sum_exp(-g[keep] %*% a[[ai]]$g[keep,] / Tfrac))[1, 1]
    
    for(j in 1:length(keep)) {
      k = keep[j]
      if(k == 14) {
        grad[[j]] = -3 * (data[[i]]$Eg[3] * data[[i]]$Eg[2] - a[[ai]]$Eg[3] * a[[ai]]$Eg[2]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(3 * a[[ai]]$g[3,] * a[[ai]]$g[2,]) * N * N
      } else if(k == 15) {
        grad[[j]] = -(data[[i]]$Eg[4] * data[[i]]$Eg[2] - a[[ai]]$Eg[4] * a[[ai]]$Eg[2]) / Tfrac -
          2 * (data[[i]]$Eg[3] * data[[i]]$Eg[2] - a[[ai]]$Eg[3] * a[[ai]]$Eg[2]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(2 * a[[ai]]$g[3,] * a[[ai]]$g[2,] + a[[ai]]$g[4,] * a[[ai]]$g[2,]) * N * N
      } else if(k == 16) {
        grad[[j]] = -3 * (data[[i]]$Eg[4] * data[[i]]$Eg[2] - a[[ai]]$Eg[4] * a[[ai]]$Eg[2]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(3 * a[[ai]]$g[4,] * a[[ai]]$g[2,]) * N * N
      } else if(k == 17) {
        grad[[j]] = -(data[[i]]$Eg[3] * data[[i]]$Eg[2] - a[[ai]]$Eg[3] * a[[ai]]$Eg[2]) / Tfrac -
          (data[[i]]$Eg[4] * data[[i]]$Eg[2] - a[[ai]]$Eg[4] * a[[ai]]$Eg[2]) / Tfrac -
          (data[[i]]$Eg[5] * data[[i]]$Eg[2] - a[[ai]]$Eg[5] * a[[ai]]$Eg[2]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(a[[ai]]$g[3,] * a[[ai]]$g[2,] +
                                              a[[ai]]$g[4,] * a[[ai]]$g[2,] +
                                              a[[ai]]$g[5,] * a[[ai]]$g[2,]) * N * N
      } else if(k == 18) {
        grad[[j]] = -3 * (data[[i]]$Eg[5] * data[[i]]$Eg[2] - a[[ai]]$Eg[5] * a[[ai]]$Eg[2]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(3 * a[[ai]]$g[5,] * a[[ai]]$g[2,]) * N * N
      } else {
        grad[[j]] = -(data[[i]]$Eg[k] - a[[ai]]$Eg[k]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(a[[ai]]$g[k,]) * N * N
      }
      
      #if(k > 10) {
      #  grad[[j]] = grad[[j]] + ifelse(g[k] > 0.0, -lambda, lambda)
      #}
    }
    
    out[[ai]]$Eg = a[[ai]]$Eg
    out[[ai]]$Egrad = rep(0, length(g))
    out[[ai]]$Egrad2 = rep(0, length(g))
    out[[ai]]$Egrad[keep] = grad
    out[[ai]]$Egrad2[keep] = grad2
  }
  
  out
}

GgradG = function(path, g, mu = NULL) {
  out = list(lp = 0, Egrad = rep(0, length(g)), Egrad2 = rep(0, length(g)))
  
  for(i in 1:length(ts)) {
    setTemperatureFraction(path, ts[i])
    a = GgradG2(path, g, data[[i]], mu)
    
    out$lp = out$lp + Reduce(`+`, map(a, ~ .$lp))
    out$Egrad = out$Egrad + Reduce(`+`, map(a, ~ .$Egrad))
    out$Egrad2 = out$Egrad2 + Reduce(`+`, map(a, ~ .$Egrad2))
  }
  
  out
}

setSupercell(path, N)
setECIs(path, ecis)
ts = seq(0.5, 5, length = 15)
Sys.time()
data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(path, ecis)
}
Sys.time()

map(data, ~ .[[8]]$Eg[c(2, 3, 14)]) %>%
  do.call(rbind, .) %>% as.tibble %>%
  mutate(corr1_2 = corr1 * corr2,
         corr1_1_1 = corr1 * corr1 * corr1) %>%
  #select(corr1_2, capprox, corr1_1_1, corr13) %>%
  mutate(t = ts) %>%
  gather(which, value, -t) %>%
  ggplot(aes(t, value)) +
  geom_point(aes(shape = which, colour = which), alpha = 0.75)

map(1:length(data), function(j) {
  map(1:length(mus), ~ data[[j]][[.]]$Eg[c(2, 3, 14)]) %>%
    do.call(rbind, .) %>%
    as.tibble %>%
    mutate(mu = mus,
           t = ts[j])
  }) %>%
  bind_rows %>%
  mutate(corr1_2 = corr1 * corr2,
         corr1_1_1 = corr1 * corr1 * corr1) %>%
  gather(which, value, -mu, -t) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(shape = which, colour = which), alpha = 0.75) +
  facet_wrap( ~ t)
