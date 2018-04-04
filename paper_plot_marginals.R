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
N = 30
keep = 2:13
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

Sys.time()
data = list()
for(i in 1:length(ts)) {
  setTemperatureFraction(path, ts[i])
  data[[i]] = runSimulation(ecis)
}
Sys.time()

UgradU = function(q) {
  qtmp = rep(0, length(ecis))
  qtmp[nonZero] = q
  a = GgradG(qtmp)
  list(u = 0.5 * t(a$Eg) %*% a$Eg,
       dudq = t(a$Eg) %*% a$Egrad[, nonZero])
}

opts3 = list()
#radford(UgradU, 1e-2, 25, ecis[nonZero])

ps = seq(-0.25, 0.25, length = 10)

ps2 = list(p1 = ps,
          p2 = ps) %>%
  as.tibble %>%
  tidyr::expand(p1, p2)

z = list()
for(i in 1:nrow(ps2)) {
  b = opts2[[1]][nonZero]
  b[[7]] = ps2[i, 1] %>% as.numeric
  b[[8]] = ps2[i, 2] %>% as.numeric
  z[[i]] = UgradU(b)
  cat("i: ", i, " out of ", nrow(ps2), "\n")
}

ps2 %>%
  mutate(z = map(z, ~ .$u)) %>%
  unnest() %>%
  mutate(z = (z < 5) * z + (z > 5) * 5) %>%
  ggplot(aes(p1, p2)) +
  geom_tile(aes(fill = z))

  mutate(z = pmap(., function(p1, p2) {
    p1
  })) %>%
  unnest()

w = getW(opts2[[1]])
a = GgradG(opts2[[1]])

opts_full2[[5]] %>%
  map(~ as.tibble(t(.))) %>%
  bind_rows %>%
  select(nonZero) %>%
  mutate(t = row_number()) %>%
  gather(corr, value, starts_with("corr")) %>%
  ggplot(aes(t, value)) +
  geom_line(aes(colour = corr))

pc = princomp(opts2 %>%
           do.call(rbind, .) %>%
           as.tibble %>%
           select(nonZero))

pc$loadings

pc$sdev

# do l1 regularization stuff
