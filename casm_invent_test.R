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

nonZero = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
keep = 3:13
makeECIs = function() {
  ecis = rnorm(length(ecis), 0.1, 0.25)
  ecis[-nonZero] = 0
  ecis
}
ecis = makeECIs()
setECIs(path, ecis)
N = 15

runMC(path)

corrs = getCorrs(path)
results = getResults(path) %>% select(param_chem_pota, everything())

corrs %>%
  group_by(mu) %>%
  filter(mci > 1) %>% ungroup() %>%
  gather(corr, value, starts_with("corr")) %>%
  ggplot(aes(mci, value)) +
  geom_line(aes(group = corr, color = mu), alpha = 0.1)


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

corrs %>% filter(mu == 5) %>% smu

runSimulation = function(g, getG = FALSE) {
  setECIs(path, g)
  runMC(path)
  corrs = getCorrs(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500) %>%
    do(out = smu(., getG)) %>% pull(out)
}

ecis = makeECIs()
Sys.time()
data = runSimulation(ecis)
Sys.time()

list(y = map(data, ~ .$Eg[[4]]) %>% unlist()) %>% as.tibble %>%
  mutate(x = row_number()) %>%
  ggplot(aes(x, y)) +
  geom_line() +
  geom_point()

# Plot curves to be fit (the data)
results = getResults(path) %>% select(param_chem_pota, everything())
results %>% select(param_chem_pota, starts_with("corr")) %>%
  gather(corr, value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, value)) +
  geom_line(aes(group = corr, colour = corr))

GgradG = function(g, getG = FALSE) {
  a = runSimulation(g, getG)

  for(i in 1:length(a)) {
    J = length(a[[i]]$Eg)
    for(j in 1:J) {
      #a[[i]]$g[j,] = a[[i]]$g[j,] - data[[i]]$Eg[j]
      a[[i]]$Eg[j] = a[[i]]$Eg[j] - data[[i]]$Eg[j]
    }
  }
  
  out = list(Eg = do.call("rbind", map(a, ~ .$Eg)),
             Egrad = do.call("rbind", map(a, ~ .$Egrad)))
  
  if(getG) {
    out$g = do.call("rbind", map(a, ~ .$g))
  }
  
  out
}

opts = list()
opts_full = list()
for(j in 1:2) {
  bs = list()
  b = rnorm(length(ecis), 0.0, 0.0) * 0
  eps = 0.1
  M = 250
  frac = 0.6
  scaling = -log(0.02) / (frac * M)
  ts = rep(0, M)
  for(i in 1:M) {
    tryCatch ({
      g = GgradG(b)
      W = diag(nrow(g$Eg))
      u = t(g$Eg) %*% solve(W, g$Eg)
      grad = (t(g$Eg) %*% solve(W, g$Egrad))[1,]
      dt = if(i < frac * M) eps * exp(-i * scaling) else dt
      ts[[i]] = dt
      
      cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, ", params : ", b, "\n")
      bs[[i]] = b
      b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }
  
  ts = cumsum(ts)
  opts[[j]] = b
  opts_full[[j]] = bs
}

ecis %>% setNames(names(opts[[1]]))
# Compare results of optimization to true ecis
do.call(rbind, opts) %>% as.tibble %>%
  mutate(which = "optimization") %>%
  bind_rows(ecis %>% setNames(names(opts[[1]]))) %>%
  mutate(which = replace(which, is.na(which), "truth")) %>%
  select(nonZero, which) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = corr %>% unique %>% mixedsort)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 3)

# Plot optimization trajectories
do.call(rbind, opts_full[[1]]) %>% as.tibble %>%
  select(nonZero) %>% mutate(t = c(ts[1:246], cumsum(rep(dt, 50)) + ts[246])) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = unique(corr))) %>%
  ggplot(aes(t, eci)) +
  geom_line(aes(group = corr, color = corr))

# Compare both optimizations
do.call(rbind, opts) %>% as.tibble %>%
  mutate(which = "optimization") %>%
  bind_rows(ecis %>% setNames(names(opts[[1]]))) %>%
  mutate(which = replace(which, is.na(which), "truth")) %>%
  select(nonZero, which) %>%
  gather(corr, eci, starts_with("corr")) %>%
  bind_rows(do.call(rbind, opts2) %>% as.tibble %>%
              mutate(which = "optimization_w_extra_ecis") %>%
              bind_rows(ecis %>% setNames(names(opts2[[1]]))) %>%
              mutate(which = replace(which, is.na(which), "truth")) %>%
              select(c(3:12, 14:23), which) %>%
              gather(corr, eci, starts_with("corr"))) %>%
  mutate(corr = factor(corr, levels = corr %>% unique %>% mixedsort)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 3)

for(j in 1:2) {
  bs = opts_full[[j]]
  b = bs[[M]]
  for(i in (M + 1):(M + 50)) {
    tryCatch ({
      g = GgradG(b)
      W = diag(nrow(g$Eg))
      u = t(g$Eg) %*% solve(W, g$Eg)
      grad = (t(g$Eg) %*% solve(W, g$Egrad))[1,]
      dt = 0.02
      
      cat("j : ", j, ", i : ", i, ", dt : ", dt, ", lp : ", u, ", params : ", b, "\n")
      b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
      bs[[i]] = b
    }, error = function(e) {
      cat("Error at ", j, " ", i, "\n")
      cat(paste(e), "\n")
    })
  }
  
  opts[[j]] = b
  opts_full[[j]] = bs
}

