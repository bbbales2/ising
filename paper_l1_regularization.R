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

path = "/home/bbales2/casm/invent2"
ecis = rep(0, length(getECIs(path)))
N = getSupercell(path)
keep = 2:13
ts = c(0.5)
nonZero = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
ecis[nonZero] = c(0.440, 0.280, 0.100, -0.100, -0.133, 0.010, -0.180, 0.170, 0.080, -0.080)

paths = c("/home/bbales2/ising/paper_outputs/test2.dat",
          "/home/bbales2/ising/paper_outputs/test3.dat")

env = new.env()
load(paths[1], envir = env)

getData = function(data, vars) {
  data = map(paths, function(path) {
    env = new.env()
    load(path, envir = env)
    out = list()
    
    for(i in 1:length(vars)) {
      if(class(env[[vars[i]]])[[1]] == "list") {
        out[[vars[i]]] = do.call(rbind, env[[vars[i]]]) %>%
          as.tibble %>%
          mutate(which = basename(path),
                 opt = row_number())
      } else {
        out[[vars[i]]] = env[[vars[i]]] %>%
          mutate(which = basename(path))
      }
    }
    
    out
  })
  
  outT = list()
  for(i in 1:length(vars)) {
    outT[[vars[i]]] = map(data, ~ .[[vars[i]]]) %>%
      do.call(rbind, .)
  }
  
  outT
}

makeECIs = function() {
  ecis[nonZero] = rnorm(length(nonZero), 0.1, 0.25)
  ecis[-nonZero] = 0
  ecis
}

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

data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(ecis)
}

UgradU = function(q) {
  qtmp = rep(0, length(ecis))
  qtmp[nonZero] = q
  a = GgradG(qtmp)
  list(u = 0.5 * t(a$Eg) %*% a$Eg,
       dudq = t(a$Eg) %*% a$Egrad[, nonZero])
}

bs = list()
b = rep(0, length(nonZero))
eps = 0.02
M = 120
frac = 1.0
scaling = -log(0.04) / (frac * M)
for(i in 1:M) {
  ugrad = UgradU(b)
  u = ugrad$u
  grad = ugrad$dudq
  dt = if(i < frac * M) eps * exp(-i * scaling) else dt
    
  cat("i : ", i, ", dt : ", dt, ", lp : ", u, "params: \n")
  print(rbind(b, ecis[nonZero]))
  bs[[i]] = b
  b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
}

ts = rep(0, M)
for(i in 1:M) {
  ts[[i]] = if(i <= frac * M) eps * exp(-i * scaling) else dt
}
plot(ts)

b_ = rep(0, length(ecis))
b_[nonZero] = b
data2 = runSimulation(b_)

tmpFunc = function(data) {
  map(data, ~ .$Eg) %>%
    do.call(cbind, .) %>%
    t %>%
    as.tibble %>%
    mutate(mu = row_number()) %>%
    gather(corr, value, starts_with("corr"))
}

bind_rows(tmpFunc(data[[1]]) %>%
            mutate(name = "data"),
          tmpFunc(data2) %>%
            mutate(name = "fit")) %>%
  bind_rows %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(color = name)) +
  facet_wrap(~ corr)
