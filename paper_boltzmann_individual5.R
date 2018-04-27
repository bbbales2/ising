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
keep = c(2, 3, 4, 5, 6, 7, 14, 15)#, 14, 15, 16, 17, 18)
ts = c(1.0)
ecis[nonZero] = c(0.0, 0.440, 0.280, 0.100, -0.100, -0.133, 0.15, -0.33)#, 0.170, 0.080, -0.080)
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
      } else {
        grad[[j]] = -(data[[i]]$Eg[k] - a[[ai]]$Eg[k]) / Tfrac
        grad2[[j]] = (-1 / (Tfrac^2)) * var(a[[ai]]$g[k,]) * N * N
      }
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

Sys.time()
data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(path, ecis)
  for(j in 1:length(data[[i]])) {
    data[[i]][[j]]$g = NULL
    data[[i]][[j]]$Eg[-keep[1:6]] = 0.0
    data[[i]][[j]]$Egrad[,-keep[1:6]] = 0.0
  }
}
Sys.time()

opts = list()
opts_full = list()
Sys.time()
for(j in 1:1) {
  bs = list()
  b = makeECIs()
  eps = 0.1
  M = 125
  frac = 0.75
  scaling = -log(0.04) / (frac * M)
  dt = eps
  print(rbind(b[nonZero], ecis[nonZero]))
  for(i in 1:M) {
    tryCatch ({
      g = GgradG(path, b)
      u = g$lp
      grad = g$Egrad
      grad2 = g$Egrad2
      
      cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, "params: \n")
      bs[[i]] = b
      b = b - grad / (grad2 + 1e-10)
      print(rbind(b[nonZero], ecis[nonZero]))
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }
  
  list(opts = b,
       opts_full = bs)
}
Sys.time()

