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

parser$add_argument("--config", default="", help = "File to be sourced that has all config stuff in it")
parser$add_argument("--output", default="", help = "Place to store the output workspace")

args = parser$parse_args()

if(length(args$config) == 0) {
    print("No configuration file provided")
    quit()
}

if(length(args$output) == 0) {
    print("No output filename provided")
    quit()
}

source(args$config)

cat(path, "\n")
cat(N, "\n")
cat(ecis, "\n")
cat(nonZero, "\n")
cat(keep, "\n")
if(exists("use_t1")) {
  ts = c()
  
  if(use_t1) {
    ts = c(ts, 1.0)
  }
  
  if(use_t2) {
    ts = c(ts, 0.5)
  }
}
cat(ts, "\n")

makeECIs = function() {
    ecis[nonZero] = rnorm(length(nonZero), 0.1, 0.25)
    ecis[-nonZero] = 0
    ecis
}

setSupercell(path, N)

smu = function(corrs, Tfrac, getG = FALSE) {
    NN = N * N
    res = corrs %>% select(starts_with("corr")) %>% as.matrix
    
    f = matrix(0, nrow = ncol(res), ncol = 1)
    rownames(f) = colnames(res)
    jac = matrix(0, nrow = ncol(res), ncol = ncol(res))
    rownames(jac) = colnames(res)
    colnames(jac) = colnames(res)
    jac = -cov(res, res) * NN / Tfrac
    for(i in 1:ncol(res)) {
        f[i, 1] = mean(res[, i])
    }
    
    jac[, -nonZero] = 0

    ## I'm not sure why I need the two transposes on Eg to get things to work, but seems like I do
    if(getG) {
        out = list(g = t(res[,keep]), Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
    } else {
        out = list(Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
    }
    
    out
}

runSimulation = function(g, getG = FALSE) {
    setECIs(path, g)
    runMC(path)
    corrs = getCorrs(path)
    
    Tfrac = getTemperatureFraction(path)
    
    corrs %>%
        group_by(mu) %>%
        filter(mci > 500) %>%
        do(out = smu(., Tfrac, getG)) %>% pull(out)
}

GgradG2 = function(g, data, getG = FALSE) {
    a = runSimulation(g, getG)
    
    for(i in 1:length(a)) {
      J = length(a[[i]]$Eg)
      for(j in 1:J) {
        a[[i]]$Eg[j] = a[[i]]$Eg[j] - data[[i]]$Eg[j]
        if(getG) {
          a[[i]]$g[j,] = a[[i]]$g[j,] - data[[i]]$Eg[j]
        }
      }
    }
    
    out = list(Eg = do.call("rbind", map(a, ~ .$Eg)),
               Egrad = do.call("rbind", map(a, ~ .$Egrad)))
    
    if(getG) {
        out$g = do.call("rbind", map(a, ~ .$g))
    }
    
    out
}

GgradG = function(g, getG = FALSE) {
    a = list()

    for(i in 1:length(ts)) {
      setTemperatureFraction(path, ts[i])
      a[[i]] = GgradG2(g, data[[i]], getG)
    }
    
    out = list(Eg = do.call("rbind", map(a, ~ .$Eg)),
                 Egrad = do.call("rbind", map(a, ~ .$Egrad)))
      
    if(getG) {
      out$g = do.call("rbind", map(a, ~ .$g))
    }

    out
}

## Fine tune the optimizations
getW = function(b) {
    g = GgradG(b, TRUE)$g
    w_inv = matrix(0, nrow = nrow(g), ncol = nrow(g))
    for(i in 1:ncol(g)) {
        w_inv = w_inv + g[, i] %*% t(g[, i])
    }
    #w_inv = w_inv + diag(nrow(g))
    solve(w_inv) / ncol(g)
}

Sys.time()
data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(ecis)
}
Sys.time()

opts = list()
opts_full = list()
for(j in 1:20) {
    bs = list()
    b = makeECIs()
    eps = 0.1
    M = 120
    frac = 1.0
    scaling = -log(0.04) / (frac * M)
    for(i in 1:M) {
        tryCatch ({
            g = GgradG(b)
            W = diag(nrow(g$Eg))
            u = t(g$Eg) %*% solve(W, g$Eg)
            grad = (t(g$Eg) %*% solve(W, g$Egrad))[1,]
            dt = if(i < frac * M) eps * exp(-i * scaling) else dt
            
            cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, "params: \n")
            print(rbind(b[nonZero], ecis[nonZero]))
            bs[[i]] = b
            b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
        }, error = function(e) {
            cat("Error at ", j, " ", i, "\n")
            cat(paste(e), "\n")
        })
    }
    
    opts[[j]] = b
    opts_full[[j]] = bs
}

save.image(args$output)

opts2 = list()
opts_full2 = list()
for(j in 1:length(opts)) {
    ##W = diag(nrow(g$Eg))
    W = getW(opts[[j]])
    b = opts[[j]]
    bs = list()
    for(i in 1:50) {
        tryCatch ({
            g = GgradG(b)
            u = t(g$Eg) %*% W %*% g$Eg
            grad = (t(g$Eg) %*% W %*% g$Egrad)[1,]
            dt = 0.02
            
            cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, "params: \n")
            print(rbind(b[nonZero], ecis[nonZero]))
            b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
            bs[[i]] = b
        }, error = function(e) {
            cat("Error at ", j, " ", i, "\n")
            cat(paste(e), "\n")
        })
    }
    
    opts2[[j]] = b
    opts_full2[[j]] = bs
}

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