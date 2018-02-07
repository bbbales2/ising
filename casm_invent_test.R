source("casm_helper.R")
source("ising_helpers.R")
require(Rcpp)

library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(stats)
library(parallel)

path = "/home/bbales2/casm/invent"
ecis = getECIs(path)
ecis = rnorm(length(ecis), 0.1, 0.25)
setECIs(path, ecis)
N = 15

runMC(path)

corrs = getCorrs(path)
results = getResults(path) %>% select(param_chem_pota, everything())

corrs %>%
  group_by(mu) %>%
  filter(mci > 500) %>% ungroup() %>%
  gather(corr, value, starts_with("corr")) %>%
  ggplot(aes(mci, value)) +
  geom_line(aes(group = corr, color = mu), alpha = 0.1)



smu = function(corrs, getG = FALSE) {
  NN = N * N
  res = corrs %>% select(starts_with("corr")) %>% as.matrix
  keep = 1:ncol(res)
  
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
  
  # I'm not sure why I need the two transposes on Eg to get things to work, but seems like I do
  if(getG) {
    out = list(g = t(res[,keep]) / NN, Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
  } else {
    out = list(Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
  }
  
  out
}

runSimulation = function(g) {
  setECIs(path, g)
  runMC(path)
  corrs = getCorrs(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500) %>%
    do(out = smu(.)) %>% pull(out)
}

ecis = rnorm(length(ecis), 0.0, 0.05)
data = runSimulation(ecis)

corrs %>% filter(mu == 1) %>% smu

list(y = map(data, ~ .$Eg[[4]]) %>% unlist()) %>% as.tibble %>%
  mutate(x = row_number()) %>%
  ggplot(aes(x, y)) +
  geom_line() +
  geom_point()

results = getResults(path) %>% select(param_chem_pota, everything())
results %>% select(param_chem_pota, starts_with("corr")) %>%
  gather(corr, value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, value)) +
  geom_line(aes(group = corr, colour = corr))

GgradG = function(g) {
  a = runSimulation(g)

  for(i in 1:length(a)) {
    J = length(a[[i]]$Eg)
    for(j in 1:J) {
      #a[[i]]$g[j,] = a[[i]]$g[j,] - data[[i]]$Eg[j]
      a[[i]]$Eg[j] = a[[i]]$Eg[j] - data[[i]]$Eg[j]
    }
  }
  
  list(#g = do.call("rbind", map(a, ~ .$g)),
       Eg = do.call("rbind", map(a, ~ .$Eg)),
       Egrad = do.call("rbind", map(a, ~ .$Egrad)))
}

opts = list()
for(j in 1:1) {
  b = rnorm(length(ecis), 0.0, 0.05) * 0
  eps = 0.1
  M = 500
  scaling = -log(0.01) / M
  for(i in 1:M) {
    g = GgradG(b)
    W = diag(nrow(g$Eg))
    u = t(g$Eg) %*% solve(W, g$Eg)
    grad = (t(g$Eg) %*% solve(W, g$Egrad))[1,]
    dt = eps * exp(-i * scaling)
    
    cat("j : ", j, ", dt : ", dt, ", lp : ", u, ", params : ", b, "\n")
    b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
  }
  
  opts[[j]] = b
}
