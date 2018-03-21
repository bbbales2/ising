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

fd(function(x) {  }, ecis[nonZero[1]], 0.1)

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
  
  Tfrac = getTemperatureFraction(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500) %>%
    do(out = smu(., Tfrac, getG)) %>% pull(out)
}

# lwrap = function(b, i) {
#   ecis2 = ecis3
#   ecis2[nonZero[i]] = b
#   setTemperatureFraction(path, 0.5)
#   runSimulation(ecis2)[[10]]
# }
# 
# dx = 1e-2
# l1 = lwrap(ecis3[nonZero[1]] + dx, 1)
# l2 = lwrap(ecis3[nonZero[1]], 1)
# l3 = lwrap(ecis3[nonZero[1]] - dx, 1)
# paste(sum(l1$Eg - l2$Eg) / dx, sum(l1$Egrad[,nonZero[1]]), sum(l2$Egrad[,nonZero[1]]))
# paste(sum(l2$Eg - l3$Eg) / dx, sum(l2$Egrad[,nonZero[1]]), sum(l3$Egrad[,nonZero[1]]))

corrs %>% filter(mu == 5) %>% smu

ecis = makeECIs()
Sys.time()
setTemperatureFraction(path, 1.0)
data1 = runSimulation(ecis)
setTemperatureFraction(path, 0.5)
data2 = runSimulation(ecis)
Sys.time()

# Plot curves to be fit (the data)
getResults(path) %>% select(param_chem_pota, everything()) %>%
  select(param_chem_pota, starts_with("corr")) %>%
  select(1, mixedsort(names(.))[keep]) %>%
  gather(corr, avg_value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, avg_value)) +
  geom_line(aes(group = corr, colour = corr)) +
  theme_grey(base_size = 18)

GgradG2 = function(g, data, getG = FALSE) {
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

GgradG = function(g, getG = FALSE) {
  setTemperatureFraction(path, 1.0)
  a1 = GgradG2(g, data1, getG)
  setTemperatureFraction(path, 0.5)
  a2 = GgradG2(g, data2, getG)
  
  out = list(Eg = rbind(a1$Eg, a2$Eg),
             Egrad = rbind(a1$Egrad, a2$Egrad))
  
  if(getG) {
    out$g = rbind(a1$g, a2$g)
  }
  
  out
}

opts = list()
opts_full = list()
for(j in 1:20) {
  bs = list()
  b = makeECIs()#rnorm(length(ecis), 0.0, 0.0) * 0#
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

# Make time
ts = rep(0, M)
for(i in 1:M) {
  dt = if(i < frac * M) eps * exp(-i * scaling) else dt
  ts[[i]] = dt
}
ts = cumsum(ts)

ecis %>% setNames(names(opts[[1]]))
# Compare results of optimization to true ecis
do.call(rbind, opts2) %>% as.tibble %>%
  mutate(which = "optimization") %>%
  bind_rows(ecis %>% setNames(names(opts3[[1]]))) %>%
  mutate(which = replace(which, is.na(which), "truth")) %>%
  select(nonZero, which) %>%
  mutate(opt = row_number()) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = corr %>% unique %>% mixedsort)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 3, stroke = 1)

# Plot optimization trajectories
i = 2
do.call(rbind, opts_full2[[i]]) %>% as.tibble %>%
  select(nonZero) %>% mutate(t = 1:50) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = unique(corr))) %>%
  ggplot(aes(t, eci)) +
  geom_line(aes(group = corr, color = corr))
rbind(opts[[i]][nonZero], ecis[nonZero])

# Compare both optimizations
do.call(rbind, opts2) %>% as.tibble %>%
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

setTemperatureFraction(path, 0.5)
i = 1

runSimulation(ecis)
datr = getResults(path) %>% mutate(data = TRUE)
runSimulation(opts2[[i]])
getResults(path) %>% mutate(data = FALSE) %>%
  bind_rows(datr) %>%
  select(data, param_chem_pota, starts_with("corr")) %>%
  select(1:2, mixedsort(names(.))[14:18]) %>%
  gather(corr, avg_value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, avg_value)) +
  geom_point(aes(group = corr, colour = corr, shape = data), alpha = 0.5, size = 2) +
  theme_grey(base_size = 18)

# Fine tune the optimizations
getW = function(b) {
  g = GgradG(b, TRUE)$g
  w_inv = matrix(0, nrow = nrow(g), ncol = nrow(g))
  for(i in 1:ncol(g)) {
    w_inv = w_inv + g[, i] %*% t(g[, i])
  }
  solve(w_inv) / ncol(g)
}

opts2 = list()
opts_full2 = list()
for(j in 1:length(opts)) {
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

opts3 = list()
opts_full3 = list()
for(j in 1:length(opts2)) {
  W = getW(opts2[[j]])
  b = opts2[[j]]
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
  
  opts3[[j]] = b
  opts_full3[[j]] = bs
}

# Make convex hull plots
tclex = getClex(path, ecis) %>% mutate(which = "truth")
clexes = list()
for(j in 1:length(opts2)) {
  clexes[[j]] = getClex(path, opts2[[j]]) %>% mutate(which = "optimization", opt = j)
}

hull = tclex %>% group_by(comp) %>% filter(row_number() == which.min(formation_energy))#top_n(-1, formation_energy)

do.call("bind_rows", clexes) %>%
  filter(configname %in% (hull %>% pull(configname))) %>%
  bind_rows(hull) %>%
  ggplot(aes(comp, formation_energy)) +
  geom_point(aes(color = which), shape = 4)

do.call("bind_rows", clexes) %>%
  group_by(comp, opt) %>%
  mutate(is_minimum = (row_number() == which.min(formation_energy))) %>%
  ungroup() %>%
  filter(configname %in% (hull %>% pull(configname))) %>%
  group_by(comp) %>%
  mutate(n = n()) %>%
  mutate(num_mins = sum(is_minimum)) %>%
  mutate(num_not_mins = n()) %>%
  ungroup() %>%
  gather(which_count, num, num_mins, num_not_mins) %>%
  ggplot(aes(comp, num)) +
  geom_bar(aes(fill = which_count), stat = "identity")

# Make the cooling run plots
tcr = coolingRun(path, ecis, seq(0.1, 1.0, length = 20)) %>%
  mutate(which = "truth")
crs = list()
for(j in 1:length(opts2)) {
  crs[[j]] = coolingRun(path, opts2[[j]], seq(0.1, 1.0, length = 20)) %>%
    mutate(which = "optimization", opt = j)
  cat("Finished cr: ", j, " of ", length(opts2), "\n")
}

crs %>% bind_rows %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_point(aes(colour = param_chem_pota), alpha = 0.1) +
  geom_path(data = tcr, aes(group = param_chem_pota), colour = "red")

crs %>% bind_rows %>%
  mutate(chem = factor(param_chem_pota, levels = sample(unique(param_chem_pota)))) %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_path(aes(colour = chem), alpha = 0.5) +
  geom_path(data = tcr, aes(group = param_chem_pota), colour = "black", alpha = 0.5)
