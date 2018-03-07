# This file is incomplete -- requires a bunch of stuff from casm_invent_test.R
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
ecis[nonZero] = c(0.2, 0.05760729, 0.13852900, 0.22172975, 0.17998039, -0.03130637, 0.18242629, 0.20327261, -0.11770949, 0.05011442)
ecis = makeECIs(); ecis
setECIs(path, ecis)
N = 45
setSupercell(path, N)

runMC(path)

corrs = getCorrs(path)
results = getResults(path) %>% select(param_chem_pota, everything())
ecis = getECIs(path)

corrs %>%
  group_by(mu) %>%
  filter(mci > 800) %>% ungroup() %>%
  gather(corr, value, corr2) %>%#starts_with("corr")
  ggplot(aes(mci, value)) +
  geom_point(aes(group = corr, color = mu), alpha = 0.1)

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
        #for(i in 1:ncol(res)) {
        #for(j in 1:ncol(res)) {
        #for(k in 1:ncol(res)) {
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

i = 10
ecis[nonZero] = c(0.2125, 0.0393, 0.1420, 0.2261, 0.1879, -0.0058, 0.1128, 0.2400, -0.0973, 0.0184)
setECIs(path, ecis)
runMC(path)
corrs = getCorrs(path)
lps = map(1:length(data), function(i) { laplace(corrs %>% filter(mci > 500, mu == i), data[[i]]$Eg) })

#do.call("sum", )
sigma = 0.01
Reduce('+', map(lps, ~ .$l / sigma^2))
Reduce('+', map(lps, ~ .$jac / sigma^2))
Reduce('+', map(lps, ~ .$hessian / sigma^2)) %>% solve -> Sigma

mvrnorm(500, ecis[nonZero], Sigma) %>%
  as.tibble %>%
  ggpairs

#lp = function(x) {
#}

mvrnorm(5, ecis[nonZero], solve(lp$hessian))

fd = function(func, x, dx = 1e-2) {
  (func(x + dx) - func(x)) / dx
}

lwrap = function(b, i) {
  ecis2 = ecis
  ecis2[nonZero[i]] = b
  setECIs(path, ecis2)
  runMC(path)
  corrs = getCorrs(path)
  laplace(corrs %>% filter(mci > 500, mu == 10), data[[10]]$Eg)
}

dx = 4e-2
l1 = lwrap(ecis2[nonZero[1]] + dx, 1)
l2 = lwrap(ecis2[nonZero[1]], 1)
l3 = lwrap(ecis2[nonZero[1]] - dx, 1)
paste((l1$l - l2$l) / dx, l1$jac[1], l2$jac[1])
paste((l1$l - 2 * l2$l + l3$l) / (dx^2), l1$hessian[1, 1], l2$hessian[1, 1])

lwrap2 = function(b1, b2, i, j) {
  ecis2 = ecis
  ecis2[nonZero[i]] = b1
  ecis2[nonZero[j]] = b2
  setECIs(path, ecis2)
  runMC(path)
  corrs = getCorrs(path)
  laplace(corrs %>% filter(mci > 500, mu == 10), data[[10]]$Eg)
}

dx = 5e-2
l11 = lwrap2(ecis2[nonZero[1]] + dx, ecis2[nonZero[2]] + dx, 1, 2)
l12 = lwrap2(ecis2[nonZero[1]] + dx, ecis2[nonZero[2]], 1, 2)
l21 = lwrap2(ecis2[nonZero[1]], ecis2[nonZero[2]] + dx, 1, 2)
l22 = lwrap2(ecis2[nonZero[1]], ecis2[nonZero[2]], 1, 2)
#paste((l1$l - l2$l) / dx, l1$jac[1], l2$jac[1])
paste(((l11$l - l12$l) - (l21$l - l22$l)) / (dx^2), l11$hessian[1, 2], l11$hessian[1, 2], l21$hessian[1, 2], l22$hessian[1, 2])
