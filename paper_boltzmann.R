#!/usr/bin/env Rscript
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

parser = ArgumentParser()

parser$add_argument("--config", default="", help = "File to be sourced that has all config stuff in it (or an output file if --refine enabled)")
parser$add_argument("--output", default="", help = "Place to store the output workspace")
parser$add_argument("--refine", action="store_true", default=FALSE, help = "Refine results instead of computing new ones (--config should be an --output file)")

args = parser$parse_args()

if(length(args$config) == 0) {
  print("No configuration file provided")
  quit()
}

if(length(args$output) == 0) {
  print("No output filename provided")
  quit()
}

if(args$refine) {
  attach(args$config)
} else {
  source(args$config)
}

cat(path, "\n")
cat(N, "\n")
cat(ecis, "\n")
cat(nonZero, "\n")
cat(keep, "\n")
cat(ts, "\n")

makeECIs = function() {
  ecis[keep] = rnorm(length(keep), 0.1, 0.25)
  ecis[-keep] = 0
  ecis
}

setSupercell(path, N)

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

runSimulation = function(g) {
  setECIs(path, g)
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

GgradG2 = function(g, data) {
  a = runSimulation(g)
  
  Tfrac = getTemperatureFraction(path)
  
  out = list()
  for(i in 1:length(a)) {
    grad = rep(0, length(keep))
    
    out[[i]] = list()
    
    out[[i]]$lp = (-data[[i]]$Eg[keep] %*% g[keep] / Tfrac - log_sum_exp(-g[keep] %*% a[[i]]$g[keep,] / Tfrac))[1, 1]
    
    for(j in 1:length(keep)) {
      k = keep[j]
      grad[[j]] = -(data[[i]]$Eg[k] - a[[i]]$Eg[k]) / Tfrac
    }
    
    out[[i]]$Eg = a[[i]]$Eg
    out[[i]]$Egrad = rep(0, length(g))
    out[[i]]$Egrad[keep] = grad
  }
  
  out
}

GgradG = function(g) {
  out = list(lp = 0, Egrad = rep(0, length(g)))
  
  for(i in 1:length(ts)) {
    setTemperatureFraction(path, ts[i])
    a = GgradG2(g, data[[i]])
    
    out$lp = out$lp + Reduce(`+`, map(a, ~ .$lp))
    out$Egrad = out$Egrad + Reduce(`+`, map(a, ~ .$Egrad))
  }
  
  out
}

Sys.time()
data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(ecis)
}
Sys.time()

if(length(nonZero) > length(keep)) {
  toPrint = nonZero
} else {
  toPrint = keep
}

opts = list()
opts_full = list()
for(j in 1:20) {
  bs = list()
  b = makeECIs()
  eps = 0.1
  M = 50
  frac = 1.0
  scaling = -log(0.04) / (frac * M)
  for(i in 1:M) {
    tryCatch ({
      g = GgradG(b)
      u = g$lp
      grad = g$Egrad
      dt = if(i < frac * M) eps * exp(-i * scaling) else dt
      
      cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, "params: \n")
      print(rbind(b[toPrint], ecis[toPrint]))
      bs[[i]] = b
      b = b + dt * grad / (sqrt(sum(grad^2)) + 1e-10)
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }

  opts[[j]] = b
  opts_full[[j]] = bs
}

opts2 = opts
opts_full2 = opts_full

save.image(args$output)

## Make convex hull data
tclex = getClex(path, ecis) %>% mutate(which = "truth")
clexes = list()
for(j in 1:length(opts2)) {
  clexes[[j]] = getClex(path, opts2[[j]]) %>%
    mutate(which = "optimization", opt = j)
}

hull = tclex %>%
  group_by(comp) %>%
  filter(row_number() == which.min(formation_energy))

## Make the cooling run data
tcr = coolingRun(path, ecis, seq(0.1, 1.0, length = 20)) %>%
  mutate(which = "truth")
crs = list()
for(j in 1:length(opts2)) {
  crs[[j]] = coolingRun(path, opts2[[j]], seq(0.1, 1.0, length = 20)) %>%
    mutate(which = "optimization", opt = j)
  cat("Finished cr: ", j, " of ", length(opts2), "\n")
}

## Generate predictions
runSimDataToDf = function(data) {
  do.call(cbind, map(data, ~ .$Eg)) %>%
    as.tibble %>%
    setNames(getChemicalPotentials(path)) %>%
    mutate(corr = as.factor(row_number())) %>%
    gather(mu, value, -corr) %>%
    mutate(mu = as.numeric(mu))
}

runAllTemps = function(ecis_) {
  a = list()
  for(i in 1:length(ts)) {
    setTemperatureFraction(path, ts[i])
    a[[i]] = runSimDataToDf(runSimulation(ecis_)) %>%
      mutate(temp = ts[i])
  }
  
  do.call("bind_rows", a) %>%
    mutate(opt = i)
}

pData = list()
for(i in 1:length(opts2)) {
  pData[[i]] = runAllTemps(opts2[[i]])
  cat("Finished ", i, " out of ", length(opts2), "\n")
}

tpData = runAllTemps(ecis) %>%
  mutate(which = "truth")

save.image(args$output)