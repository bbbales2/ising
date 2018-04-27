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

parser = ArgumentParser()

parser$add_argument("--output", default="", help = "Place to store the output workspace")

args = parser$parse_args()

if(length(args$output) == 0) {
  print("No output filename provided")
  quit()
}

path = "/home/bbales2/casm/invent3"
N = 30
ecis = rep(0, length(getECIs(path)))
nonZero = c(2, 3, 4, 5, 6, 7, 14, 15)
keep = c(2, 3, 4, 5, 6, 7)
ts = c(4.0)
ecis[nonZero[1]] = 0.0
ecis[nonZero[2]] = rnorm(1, 0.30, 0.1)
ecis[nonZero[3:length(nonZero)]] = rnorm(length(nonZero) - 2, 0.1, 0.25)
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
    
    if(!is.null(mu)) {
      ai = 1
    } else {
      ai = i
    }
    
    out[[ai]] = list()
    
    out[[ai]]$lp = (-data[[i]]$Eg[keep] %*% g[keep] / Tfrac - log_sum_exp(-g[keep] %*% a[[ai]]$g[keep,] / Tfrac))[1, 1]
    
    for(j in 1:length(keep)) {
      k = keep[j]
      grad[[j]] = -(data[[i]]$Eg[k] - a[[ai]]$Eg[k]) / Tfrac
    }
    
    out[[ai]]$Eg = a[[ai]]$Eg
    out[[ai]]$Egrad = rep(0, length(g))
    out[[ai]]$Egrad[keep] = grad
  }
  
  out
}

GgradG = function(path, g, mu = NULL) {
  out = list(lp = 0, Egrad = rep(0, length(g)))
  
  for(i in 1:length(ts)) {
    setTemperatureFraction(path, ts[i])
    a = GgradG2(path, g, data[[i]], mu)
    
    out$lp = out$lp + Reduce(`+`, map(a, ~ .$lp))
    out$Egrad = out$Egrad + Reduce(`+`, map(a, ~ .$Egrad))
  }
  
  out
}

Sys.time()
data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(path, ecis)
}
Sys.time()

opts = list()
opts_full = list()
Sys.time()
res = mclapply(1:length(mus), function(j) {
  rpath = paste0("/home/bbales2/casm/invent", 2 + j)
  bs = list()
  b = makeECIs()
  eps = 0.1
  M = 125
  frac = 0.75
  scaling = -log(0.04) / (frac * M)
  dt = eps
  for(i in 1:M) {
    tryCatch ({
      g = GgradG(rpath, b, mus[j])
      u = g$lp
      grad = g$Egrad
      dt = if(i < frac * M) eps * exp(-i * scaling) else dt
      
      cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, "params: \n")
      print(rbind(b[nonZero], ecis[nonZero]))
      bs[[i]] = b
      b = b + dt * grad / (sqrt(sum(grad^2)) + 1e-10)
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }
  
  list(opts = b,
       opts_full = bs)
}, mc.cores = 20)
Sys.time()

opts = map(res, ~.$opts)
opts_full = map(res, ~.$opts_full)

b2 = map(opts, ~ .[[keep[2]]]) %>% unlist
b3 = map(opts, ~ .[[keep[3]]]) %>% unlist

y2 = map(data[[1]], ~ .$Eg[[2]]) %>% unlist

fit = stan("models/varying_coefficients.stan",
           data = list(N = length(y2),
                       b1 = y2,
                       b2 = b2,
                       b3 = b3),
           iter = 2000, chains = 4, cores = 4)

save.image(args$output)