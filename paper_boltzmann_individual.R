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

#there are n different pairs from which to build the triplet, so three opportunities
path = "/home/bbales2/casm/invent3"
N = 30
ecis = rep(0, length(getECIs(path)))
nonZero = c(2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
keep = c(2, 3, 4, 5, 6, 7)#, 14, 15, 16, 17, 18)
ts = c(4.0)
ecis[nonZero] = c(0.0, 0.440, 0.280, 0.100, -0.100, -0.133, 0.40, -0.1, 0, 0, 0)# 0.170, 0.080, -0.080)
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
#for(j in 1:length(mus)) {
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

opts %>%
  map(~ .[keep]) %>%
  do.call(rbind, .) -> b

data[[1]] %>%
  map(~ .$Eg[keep]) %>%
  do.call(rbind, .) -> phi

fit = stan("models/varying_coefficients_everything.stan",
           data = list(N = nrow(b),
                       M = ncol(b),
                       phi = phi,
                       b = b),
           iter = 2000, chains = 4, cores = 4)

e1 = ecis[nonZero[-1]]
e2 = extract(fit, pars = c('v', 'v13', 'v14', 'v15', 'v16', 'v17')) %>%
  do.call(cbind, .) %>%
  `[`(,-1) %>%
  colMeans %>%
  as.vector

print(rbind(e1, e2))

b %>%
  as.tibble %>%
  mutate(mu = mus) %>%
  gather(key, value, -mu) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = key))

phi %>%
  as.tibble %>%
  mutate(mu = mus) %>%
  gather(key, value, -mu) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = key))

plot(y2, b2)
# 0.4910892 0.2009558 0.1440032 -0.09527624 -0.141923
fit = stan("models/varying_coefficients.stan",
           data = list(N = length(y2),
                       b1 = y2,
                       b2 = b2,
                       b3 = b3),
           iter = 2000, chains = 4, cores = 4)

plot(map(data[[1]], ~ .$Eg[[nonZero[1]]]) %>% unlist)

b = makeECIs()
for(i in 1:length(nonZero)) {
  ind = nonZero[i]
  dx = 0.025
  b2 = b
  b2[ind] = b[ind] + dx
  l1 = GgradG(b)
  l2 = GgradG(b2)
  cat("approx grad: ", (l2$lp - l1$lp) / dx, ", grad 1: ", l1$Egrad[ind], ", grad 2: ", l2$Egrad[ind], "\n")
}

eps = 0.1
M = 125
frac = 0.75
scaling = -log(0.04) / (frac * M)
ts2 = rep(0, M)
for(i in 1:M) {
  ts2[[i]] = if(i < frac * M) eps * exp(-i * scaling) else dt
}
ts2 = cumsum(ts2)

opts_full[[2]] %>%
  do.call(rbind, .) %>%
  as.tibble %>%
  setNames(names(data[[1]][[1]]$Eg)) %>%
  select(keep) %>%
  mutate(t = ts2) %>%
  gather(corr, value, starts_with("corr")) %>%
  ggplot(aes(t, value)) +
  geom_line(aes(colour = corr, group = corr))

a = runSimulation(b)

runSimDataToDf = function(data) {
  #do.call(cbind, map(data, ~ .$Eg[c(2, 3, 14)])) %>%
  #  rbind(., .[1,] * .[2,]) %>%
  #  rbind(., .[3,] - .[4,]) %>%
  do.call(cbind, map(data, ~ .$Eg[c(2, 3, 14)])) %>%
    as.tibble %>%
    setNames(getChemicalPotentials(path)) %>%
    mutate(corr = as.factor(row_number())) %>%
    gather(mu, value, -corr) %>%
    mutate(mu = as.numeric(mu))
}

runSimDataToDf2 = function(data) {
  #do.call(cbind, map(data, ~ .$Eg[c(2, 3, 14)])) %>%
  #  rbind(., .[1,] * .[2,]) %>%
  #  rbind(., .[3,] - .[4,]) %>%
  do.call(cbind, map(data, ~ .$Ega[c(2, 3, 14)])) %>%
    as.tibble %>%
    setNames(getChemicalPotentials(path)) %>%
    mutate(corr = as.factor(row_number())) %>%
    gather(mu, value, -corr) %>%
    mutate(mu = as.numeric(mu))
}

list(runSimDataToDf(a) %>%
       mutate(which = "approx"),
     runSimDataToDf2(data[[1]]) %>%
       mutate(which = "true")) %>%
  bind_rows %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = which)) +
  facet_grid(. ~ corr)

for(j in 1:length(data[[1]])) {
  data[[1]][[j]]$Ega = data_orig[[1]][[j]]$Eg
  data[[1]][[j]]$Ega[-keep] = 0.0
  data[[1]][[j]]$Ega[14] = a[[j]]$Eg[14]
  data[[1]][[j]]$Ega[15] = a[[j]]$Eg[15]
  data[[1]][[j]]$Ega[16] = a[[j]]$Eg[16]
  data[[1]][[j]]$Ega[17] = a[[j]]$Eg[17]
  data[[1]][[j]]$Ega[18] = a[[j]]$Eg[18]
}

do.call(cbind, map(a, ~ .$Eg[c(2, 3, 14)])) %>%
  as.tibble

setECIs(path, ecis)
runMC(path)
corrs = getCorrs(path)

Tfrac = getTemperatureFraction(path)

corrs %>%
  group_by(mu) %>%
  filter(mci > 500) %>%
  select(corr1, corr3, corr14) %>%
  summarize(actual = mean(corr14),
            triplet = mean(corr1 * corr3),
            first = mean(corr1) * mean(corr3),
            corr1 = mean(corr1),
            corr2 = mean(corr3)) %>%
  gather(moment, value, -mu) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(color = moment, shape = moment))
