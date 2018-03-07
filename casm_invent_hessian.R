library(MASS)
source("casm_helper.R")
source("ising_helpers.R")
require(Rcpp)

library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(stats)
library(parallel)
library(gtools)
library(GGally)

path = "/home/bbales2/casm/invent"
ecis = getECIs(path)

nonZero = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
keep = 2:13
makeECIs = function() {
  ecis[nonZero] = rnorm(length(ecis), 0.1, 0.25)#c(rnorm(1, 0.1, 0.25), rnorm(4, 0.05, 0.15), rnorm(5, 0.05, 0.1))
  ecis[-nonZero] = 0
  ecis
}

N = 30
setSupercell(path, N)

ecis[nonZero] = c(0.2, 0.05760729, 0.13852900, 0.22172975, 0.17998039, -0.03130637, 0.18242629, 0.20327261, -0.11770949, 0.05011442)
ecis = makeECIs(); ecis
setECIs(path, ecis)

smu = function(corrs, getG = FALSE) {
  NN = N * N
  res = corrs %>% select(starts_with("corr")) %>% as.matrix
  
  f = matrix(0, nrow = ncol(res), ncol = 1)
  rownames(f) = colnames(res)
  jac = matrix(0, nrow = ncol(res), ncol = ncol(res))
  rownames(jac) = colnames(res)
  colnames(jac) = colnames(res)
  jac = -cov(res, res) * NN
  for(i in 1:ncol(res)) {
    f[i, 1] = mean(res[, i])
    #for(j in 1:ncol(res)) {
    #  jac[i, j] = -cov(res[, i], res[, j]) * NN
    #}
  }
  
  jac[, -nonZero] = 0
  
  # I'm not sure why I need the two transposes on Eg to get things to work, but seems like I do
  if(getG) {
    out = list(g = t(res[,keep]) / NN, Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
  } else {
    out = list(Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
  }
  
  out
}
runSimulation = function(g, getG = FALSE) {
  setECIs(path, g)
  runMC(path)
  corrs = getCorrs(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500) %>%
    do(out = smu(., getG)) %>% pull(out)
}

data = runSimulation(ecis)

# Plot curves to be fit (the data)
getResults(path) %>% select(param_chem_pota, everything()) %>%
  select(param_chem_pota, starts_with("corr")) %>%
  select(1, mixedsort(names(.))[keep]) %>%
  gather(corr, avg_value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, avg_value)) +
  geom_line(aes(group = corr, colour = corr)) +
  theme_grey(base_size = 18)

laplace = function(corrs, data) {
  NN = N * N
  res = corrs %>% dplyr::select(starts_with("corr")) %>% as.matrix
  res = res * NN
  
  f = array(0, c(length(keep)))
  rownames(f) = colnames(res[keep])
  jac = array(0, c(length(keep), length(nonZero)))
  rownames(jac) = colnames(res[keep])
  colnames(jac) = colnames(res[keep])
  #jac = -cov(res, res) * NN
  hessian = array(0, c(length(keep), length(nonZero), length(nonZero)))
  for(ii in 1:length(keep)) {
    i = keep[ii]
    f[ii] = mean(res[, i]) / NN
    for(jj in 1:length(nonZero)) {
      j = nonZero[jj]
      jac[ii, jj] = -cov(res[, i], res[, j]) / NN
      for(kk in 1:length(nonZero)) {
        k = nonZero[kk]
        hessian[ii, jj, kk] = (mean(res[, i] * res[, j] * res[, k]) -
                                 mean(res[, i] * res[, j]) * mean(res[, k]) -
                                 cov(res[, i], res[, k]) * mean(res[, j]) -
                                 cov(res[, j], res[, k]) * mean(res[, i])) / NN
      }
    }
  }
  
  l = 0.0
  ljac = array(0, c(1, length(nonZero)))
  lhessian = array(0, c(length(nonZero), length(nonZero)))
  for(i in 1:length(data)) {
    l = l + (f[i] - data[[i]])^2
    for(j in 1:length(nonZero)) {
      ljac[j] = ljac[j] + 2 * (f[i] - data[[i]]) * jac[i, j]
      for(k in 1:length(nonZero)) {
        lhessian[j, k] = lhessian[j, k] + 2 * jac[i, j] * jac[i, k] + 2 * (f[i] - data[[i]]) * hessian[i, j, k]
      }
    }
  }
  
  list(l = l, jac = ljac, hessian = lhessian)
}

getLp = function(b) {
  setECIs(path, b)
  runMC(path)
  corrs = getCorrs(path)
  lps = map(1:length(data), function(i) { laplace(corrs %>% filter(mci > 500, mu == i), data[[i]]$Eg) })

  sigma = 0.01
  list(l = Reduce('+', map(lps, ~ .$l / sigma^2)),
       jac = Reduce('+', map(lps, ~ .$jac / sigma^2)),
       hessian = Reduce('+', map(lps, ~ .$hessian / sigma^2)))
}

for(j in 1:1) {
  bs = list()
  b = rnorm(length(ecis), 0.0, 0.0) * 0#makeECIs()#ecis#
  eps = 0.1
  M = 120
  frac = 1.0
  scaling = -log(0.04) / (frac * M)
  for(i in 1:M) {
    tryCatch ({
      lp = getLp(b)
      dt = if(i < frac * M) eps * exp(-i * scaling) else dt

      u = lp$l
      grad = lp$jac
      cat("j : ", j, ", i : ", i, ", lp : ", u, "params: \n")
      print(rbind(b[nonZero], ecis[nonZero]))
      bs[[i]] = b
      b[nonZero] = b[nonZero] - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }
  
  opts[[j]] = b
  opts_full[[j]] = bs
}

b1 = b

for(j in 1:1) {
  for(i in 1:M) {
    tryCatch ({
      lp = getLp(b1)
      
      u = lp$l
      grad = lp$jac
      cat("j : ", j, ", i : ", i, ", lp : ", u, "params: \n")
      print(rbind(b1[nonZero], ecis[nonZero]))
      bs[[i]] = b1
      b1[nonZero] = b1[nonZero] - solve(lp$hessian, as.vector(lp$jac))
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }
}

mvrnorm(500, ecis[nonZero], Sigma) %>%
  as.tibble %>%
  ggpairs
