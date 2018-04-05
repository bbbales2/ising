library(MASS)

source("casm_helper.R")
source("ising_helpers.R")
require(Rcpp)

library(reshape2)
library(tidyverse)
library(stats)
library(parallel)
library(gtools)
library(ggplot2)
library(rstan)
library(shinystan)

path = "/home/bbales2/casm/invent2"
ecis = rep(0, length(getECIs(path)))
N = getSupercell(path)
keep = 2:13
ts = c(1.0)
nonZero = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
ecis[nonZero] = c(0.440, 0.280, 0.100, -0.100, -0.133, 0.010, -0.180, 0.170, 0.080, -0.080)
#nonZero = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
#c(0.200, 0.058, 0.139, 0.222, 0.180, -0.031, 0.182, 0.203, -0.118, 0.050)
paths = c("/home/bbales2/ising/paper_outputs/test4.dat")

env = new.env()
load(paths[1], envir = env)

smu = function(corrs, Tfrac, getG = FALSE) {
  NN = N * N
  res = corrs %>% dplyr::select(starts_with("corr")) %>% as.matrix
  res = res * NN
  
  f = array(0, c(length(keep)))
  rownames(f) = colnames(res[keep])
  jac = array(0, c(length(keep), length(nonZero)))
  rownames(jac) = colnames(res[keep])
  colnames(jac) = colnames(res[keep])
  #jac = -cov(res, res) * NN
  for(ii in 1:length(keep)) {
    i = keep[ii]
    f[ii] = mean(res[, i]) / NN
    for(jj in 1:length(nonZero)) {
      j = nonZero[jj]
      jac[ii, jj] = -cov(res[, i], res[, j]) / (Tfrac * NN)
    }
  }
  
  out = list(Eg = f, Egrad = jac)
  
  if(getG) {
    out$g = t(res) / NN
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

data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(ecis)
}

cdata = map(do.call(c, data), ~ .$Eg) %>%
  do.call(cbind, .)

model = stan_model("models/linear_expansion.stan")

iterateECIs = function(ecis_old) {
  ecis_new = list()
  fits = list()
  for(i in 1:length(ecis_old)) {
    seci = rep(0, length(ecis))
    seci[nonZero] = ecis_old[[i]]

    corrs = list()
    for(j in 1:length(ts)) {
      setTemperatureFraction(path, ts[j])
      corrs[[j]] = runSimulation(seci)#ecis
    }
    corrs = do.call(c, corrs)

    d = list(M = length(corrs), # Number of chemical potentials
             P = length(nonZero), # Number of parameters
             N = dim(corrs[[1]]$g)[2], # Number of samples in approx. distribution
             Q = length(keep), # Number of measured moments
             corrDataEgrad = aperm(simplify2array(map(corrs, ~ t(.$Egrad))), c(3, 2, 1)),
             corrDataEg = aperm(simplify2array(map(corrs, ~ t(.$Eg))), c(3, 2, 1))[,,1],
             cdata = t(cdata),
             g0 = seci[nonZero])

    fits[[i]] = sampling(model, data = d, cores = 4)
    opt = optimizing(model, data = d)
    cat("linear i: ", i, " out of ", length(ecis_old), "\n")
    print(rbind(opt$par, ecis[nonZero], seci[nonZero]))
    ecis_new[[i]] = opt$par[1:length(nonZero)]
  }

  list(ecis = ecis_new,
       fits = fits)
}

initialECIs = map(env$opts2, ~ .[nonZero])#map(env$opts_full, ~ .[[1]])
opts3 = iterateECIs(initialECIs)
opts4 = iterateECIs(opts3)
opts5 = iterateECIs(opts4)

levels = names(env$opts[[1]][nonZero])
bind_rows(initialECIs %>%
            do.call(rbind, .) %>%
            as.tibble %>%
            mutate(type = "abc_original_optimization"),
          opts3 %>%
            do.call(rbind, .) %>%
            as.tibble %>%
            setNames(levels) %>%
            mutate(type = "refinement_first"),
          ecis[nonZero] %>%
            setNames(levels) %>%
            t %>%
            as.tibble %>%
            mutate(type = "TRUTH")) %>%
  gather(corr, value, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = levels)) %>%
  ggplot(aes(corr, value)) +
  geom_point(aes(colour = type), shape = 4, size = 2, stroke = 2, position = position_dodge(width = 0.75))

launch_shinystan(fit)
