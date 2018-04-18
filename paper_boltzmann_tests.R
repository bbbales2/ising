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

path = "/home/bbales2/casm/invent"
N = 30
ecis = rep(0, length(getECIs(path)))
nonZero = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
keep = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
ts = c(0.5)
ecis[nonZero] = c(0.440, 0.280, 0.100, -0.100, -0.133, 0.010, -0.180, 0.170, 0.080, -0.080)

makeECIs = function() {
  ecis[nonZero] = rnorm(length(nonZero), 0.1, 0.25)
  ecis[-nonZero] = 0
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
      grad[[j]] = -(data[[i]]$Ega[k] - a[[i]]$Eg[k]) / Tfrac
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
mus = getChemicalPotentials(path)
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(ecis)
  for(j in 1:length(data[[i]])) {
    if(j == 1) {
      data[[i]][[j]]$dEgdu = (data[[i]][[j + 1]]$Eg - data[[i]][[j]]$Eg) / (mus[j + 1] - mus[j])
    } else if(j == length(data[[i]])) {
      data[[i]][[j]]$dEgdu = (data[[i]][[j]]$Eg - data[[i]][[j - 1]]$Eg) / (mus[j] - mus[j - 1])
    } else {
      data[[i]][[j]]$dEgdu = (data[[i]][[j + 1]]$Eg - data[[i]][[j - 1]]$Eg) / (mus[j + 1] - mus[j - 1])
    }
    data[[i]][[j]]$Ega = data[[i]][[j]]$Eg
    data[[i]][[j]]$Ega[-keep] = 0.0
    data[[i]][[j]]$Ega[14] = data[[i]][[j]]$Eg[2] * data[[i]][[j]]$Eg[3] - data[[i]][[j]]$dEgdu[3] / (2 * ts[i])
    data[[i]][[j]]$Ega[15] = data[[i]][[j]]$Eg[2] * data[[i]][[j]]$Eg[4] - data[[i]][[j]]$dEgdu[4] / (2 * ts[i])
    data[[i]][[j]]$Ega[16] = data[[i]][[j]]$Eg[2] * data[[i]][[j]]$Eg[4] - data[[i]][[j]]$dEgdu[4] / (2 * ts[i])
    data[[i]][[j]]$Ega[17] = data[[i]][[j]]$Eg[2] * data[[i]][[j]]$Eg[5] - data[[i]][[j]]$dEgdu[5] / (2 * ts[i])
    data[[i]][[j]]$Ega[18] = data[[i]][[j]]$Eg[2] * data[[i]][[j]]$Eg[5] - data[[i]][[j]]$dEgdu[5] / (2 * ts[i])
  }
}
Sys.time()

rbind(map(data[[1]], ~ .$Ega[14]) %>% unlist,
      map(data[[1]], ~ .$dEgdu[3]) %>% unlist)

list(runSimDataToDf2(data[[1]]) %>%
       mutate(which = "approx"),
     runSimDataToDf(data[[1]]) %>%
       mutate(which = "true")) %>%
  bind_rows %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = which)) +
  facet_grid(. ~ corr)

opts = list()
opts_full = list()
for(j in 1:1) {
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
      print(rbind(b[nonZero], ecis[nonZero]))
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

ts = rep(0, M)
for(i in 1:M) {
  ts[[i]] = if(i < frac * M) eps * exp(-i * scaling) else dt
}
ts = cumsum(ts)

opts_full[[1]] %>%
  do.call(rbind, .) %>%
  as.tibble %>%
  setNames(names(data[[1]][[1]]$Eg)) %>%
  select(nonZero) %>%
  mutate(t = ts) %>%
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
