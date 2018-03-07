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

fd(function(x) {  }, ecis[nonZero[1]], 0.1)

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

corrs %>% filter(mu == 5) %>% smu

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
  select(1, mixedsort(names(.))[keep]) %>%
  gather(corr, avg_value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, avg_value)) +
  geom_line(aes(group = corr, colour = corr)) +
  theme_grey(base_size = 18)

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
  bind_rows(ecis %>% setNames(names(opts[[1]]))) %>%
  mutate(which = replace(which, is.na(which), "truth")) %>%
  select(nonZero, which) %>%
  mutate(opt = row_number()) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = corr %>% unique %>% mixedsort)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 3, stroke = 1)

# Plot optimization trajectories
i = 1
do.call(rbind, opts_full[[i]]) %>% as.tibble %>%
  select(nonZero) %>% mutate(t = ts) %>%
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
  crs[[j]] = coolingRun(path, opts2[[j]], seq(0.1, 1.0, length = 50)) %>%
    mutate(which = "optimization", opt = j)
  cat("Finished cr: ", j, " of ", length(opts2), "\n")
}

crs %>% bind_rows %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_point(aes(colour = param_chem_pota), alpha = 0.1) +
  geom_path(data = tcr, aes(group = param_chem_pota), colour = "red")

#
i = 1
setECIs(path, opts2[[i]])
runMC(path)
corrs = getCorrs(path)
lps = map(1:length(data), function(i) { laplace(corrs %>% filter(mci > 500, mu == i), data[[i]]$Eg) })

sigma = 0.01
Reduce('+', map(lps, ~ .$l / sigma^2))
Reduce('+', map(lps, ~ .$jac / sigma^2))
Reduce('+', map(lps, ~ .$hessian / sigma^2)) %>% solve -> Sigma

point_plots <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_point(..., alpha = 0.1)
}

mvrnorm(500, opts2[[i]][nonZero], Sigma) %>%
  as.tibble %>%
  ggpairs(lower = list(continuous = point_plots))
