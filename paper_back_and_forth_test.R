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
    out = list(g = t(res), Eg = t(t(f[keep, 1])), Egrad = jac)
  } else {
    out = list(Eg = t(t(f[keep, 1])), Egrad = jac)
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

cdata = map(data[[1]], ~ .$Eg) %>%
  do.call(cbind, .)

sampleCorrs = function(g) {
  setECIs(path, g)
  runMC(path)
  corrs = getCorrs(path)
  
  Tfrac = getTemperatureFraction(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500)
}

# Stolen from https://gist.github.com/aufrank/83572
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

corrs = runSimulation(opts2[[2]], getG = TRUE)#ecis

corrs = sampleCorrs(ecis)
corrsProb = corrs %>%
  select(starts_with("corr"), mu) %>%
  select(nonZero, mu)
corrsData = corrs %>%
  select(starts_with("corr"), mu) %>%
  select(keep, mu)
#corrs2 = corrs %>%
#  select(mu, rownames(cdata), -corr1)

mus = corrs$mu %>% unique

corrProbMats = list()
corrDataMats = list()
for(i in mus) {
  corrProbMats[[i]] = corrsProb %>%
    filter(mu == mus[[i]]) %>%
    ungroup() %>%
    select(-mu) %>%
    as.matrix
  
  corrDataMats[[i]] = corrsData %>%
    filter(mu == mus[[i]]) %>%
    ungroup() %>%
    select(-mu) %>%
    as.matrix
}

loss = function(g) {
  l = 0.0
  Tfrac = getTemperatureFraction(path)
  g2 = ecis[nonZero]
  g2[1:10] = g

  for(i in mus) {
    p = softmax(-g2 %*% corrs[[i]]$g[nonZero,] / Tfrac)
    l = l + 0.5 * sum((p %*% t(corrs[[i]]$g[keep,]) + (g2 - ecis[nonZero]) %*% t(corrs[[i]]$Egrad[,nonZero]) - cdata[,i])^2)
    #rbind(t(p) %*% corrDataMats[[i]] - cdata[,i]))
  }

  l
}

x = seq(-0.5, 0.5, length = 100)
plot(x, map(x, loss))

opt = optim(rep(0, 10), loss, hessian = TRUE)#
rbind(opt$par,
      ecis[nonZero])


ecis2 = rep(0, length(ecis))
ecis2[nonZero] = opt$par
test = runSimulation(ecis2)

tdata = map(test, ~ .$Eg) %>%
  do.call(cbind, .)

tmpFunc = function(data) {
  data %>%
    t %>%
    as.tibble %>%
    mutate(mu = row_number()) %>%
    gather(corr, value, starts_with("corr"))
}

bind_rows(tmpFunc(cdata) %>%
            mutate(name = "data"),
          tmpFunc(tdata) %>%
            mutate(name = "fit")) %>%
  bind_rows %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(color = name)) +
  facet_wrap(~ corr)

model = stan_model("models/surrogate_model.stan")

d = list(M = length(corrs), # Number of chemical potentials
         P = length(nonZero), # Number of parameters
         N = dim(corrs[[1]]$g)[2], # Number of samples in approx. distribution
         Q = length(keep), # Number of measured moments 
         Tfrac = getTemperatureFraction(path),
         corrProb = aperm(simplify2array(map(corrs, ~ t(.$g[nonZero,]))), c(3, 1, 2)),
         corrData = aperm(simplify2array(map(corrs, ~ t(.$g[keep,]))), c(3, 2, 1)),
         corrDataEgrad = aperm(simplify2array(map(corrs, ~ t(.$Egrad[keep,nonZero]))), c(3, 2, 1)),
         corrDataEg = aperm(simplify2array(map(corrs, ~ t(.$Eg))), c(3, 2, 1))[,,1],
         cdata = t(cdata),
         g0 = opts2[[2]][nonZero])

fit = sampling(model,
               data = d,#ecis
         cores = 1, chains = 1, iter = 2000)

launch_shinystan(fit)

matrix[N, P] corrProb[M];
matrix[Q, N] corrData[M];
matrix[Q, P] corrEgrad[M];
vector[Q] cdata[M];
vector[P] g0;