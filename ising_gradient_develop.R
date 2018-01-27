library(argparse)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(parallel)
require(Rcpp)
library(GGally)
require(gtools)

beta = c(0.00484411, 0.26854678, -0.22950556, 0.26727645, 0.28048656)
gamma = c(-0.1238928)
names(beta) = c("b0", "b1", "b2", "b3", "b4")
names(gamma) = c("b5")

useFuncLoss = TRUE
useCorrLoss = 5

sourceCpp("covariance2.cpp")
source("ising_helpers3.R")

p = function(x) {
  as.tibble(melt(x)) %>%
    ggplot(aes(x = Var1, y = Var2, fill = factor(value))) +  
    geom_tile() +
    xlab("x") +
    ylab("y") +
    coord_equal() +
    scale_y_reverse()
}

N = 2 * 3 * 5
sigma = 0.01
S = 200
mus = seq(-5.0, 5.0, length = 21)

x = matrix(sample(c(-1, 1), N * N, replace = TRUE), nrow = N)

smu = function(mu, beta, gamma) {
  keep = 1:(1 + useCorrLoss)
  res = ising_gibbs(x, mu, beta, gamma, S, seed)$states[(S / 2 + 1) : S, ]
  
  f = matrix(0, nrow = ncol(res), ncol = 1)
  rownames(f) = colnames(res)
  jac = matrix(0, nrow = ncol(res), ncol = ncol(res) - 1)
  rownames(jac) = colnames(res)
  colnames(jac) = c(names(beta), names(gamma))
  for(i in 1:ncol(res)) {
    f[i, 1] = mean(res[, i]) / NN
    for(j in 2:ncol(res)) {
      jac[i, j - 1] = -cov(res[, i], res[, j]) / NN
    }
  }
  
  # I'm not sure why I need the two transposes on Eg to get things to work, but seems like I do
  list(g = t(res[,keep]) / NN, Eg = t(t(f[keep, 1])), Egrad = jac[keep,])
}

#GgradG = function(g) {
  beta = g[1:5]
  gamma = g[6]

  a = mclapply(mus[1:2], function(mu) { smu(mu, beta, gamma) }, mc.cores = min(length(mus), 24))

  for(i in 1:length(data)) {
    S = length(data[[i]]$f)
    for(j in 1:S) {
      a[[i]]$g[j,] = a[[i]]$g[j,] - data[[i]]$f[j]
      a[[i]]$Eg[j] = a[[i]]$Eg[j] - data[[i]]$f[j]
    }
  }

  list(g = do.call("rbind", map(a, ~ .$g)),
       Eg = do.call("rbind", map(a, ~ .$Eg)),
       Egrad = do.call("rbind", map(a, ~ .$Egrad)))
#}

#%>% do.call("rbind", .)

ising_sweep = function(x, q, S, seeds) {
  list(mu = mus,
       seed = seeds) %>%
    expand.grid %>%
    as.tibble %>%
    (function(df) split(df, 1:nrow(df))) %>%
    mclapply(function(row) {
      beta = q[1:5]
      gamma = q[6]
      res = ising_gibbs_derivs(x, row$mu, beta, gamma, S, row$seed)
      
      res$mu = row$mu
      res$seed = row$seed
      
      res
    }, mc.cores = 24)
}

data = map(ising_sweep(x, c(beta, gamma), S * 10, 2), ~ list(f = .$f))

UgradU = function(q) {
  u = 0.0
  dudq = rep(0, length(q))
  r = sample(10000000, 1)
  
  calc = ising_sweep(x, q, S, r)
  
  for(i in 1:length(calc)) {
    if(useFuncLoss) {
      u = u + 0.5 * (calc[[i]]$f[1] - data[[i]]$f[1])^2
      dudq = dudq + (calc[[i]]$f[1] - data[[i]]$f[1]) * calc[[i]]$jac[1,]
    }
    
    if(useCorrLoss > 0) {
      u = u + 0.5 * sum((calc[[i]]$f[2:(1 + useCorrLoss)] - data[[i]]$f[2:(1 + useCorrLoss)])^2)
      dudq = dudq + (calc[[i]]$f[2:(1 + useCorrLoss)] - data[[i]]$f[2:(1 + useCorrLoss)]) %*% calc[[i]]$jac[2:(1 + useCorrLoss),]
    }
  }
  
  list(u = 100 * u,
       dudq = 100 * dudq)
}

GgradG = function(q) {
  r = sample(10000000, 1)
  
  calc = ising_sweep(x, q, S, r)
  
  g = matrix(0, nrow = length(calc) * (1 + useCorrLoss), ncol = 1)
  dgdq = matrix(0, nrow = length(calc) * (1 + useCorrLoss), ncol = length(q))
  ii = 1
  for(i in 1:length(calc)) {
    for(j in 1:(1 + useCorrLoss)) {
      g[ii] = 0.5 * (calc[[i]]$f[j] - data[[i]]$f[j])^2
      dgdq[ii,] = (calc[[i]]$f[j] - data[[i]]$f[j]) * calc[[i]]$jac[j,]
      ii = ii + 1
    }
  }
  
  list(g = g,
       dgdq = dgdq)
}

fn = function(b) { UgradU(b)$u }
gn = function(b) { UgradU(b)$dudq }

b = rnorm(6, 0.1, 0.25)
eps = 10.0
for(i in 1:100) {
  ugradu = UgradU(b)
  grad = ugradu$dudq
  
  cat("lp : ", ugradu$u, ", params : ", b, "\n")
  b = b - eps * grad / (sum(grad^2) + 1e-10)
}
cat("lp : ", ugradu$u, ", params : ", b, "\n")

opts = list()
for(o in 1:400) {
  b = rnorm(6, 0.1, 0.25)
  names(b) = c("b0", "b1", "b2", "b3", "b4", "b5")
  out = optim(b, fn, gn, method = "L-BFGS-B", control = list(maxit = 200, trace = 1, REPORT = 1))
  opts[[o]] = out$par %>% as.list %>% as.tibble %>%
    mutate(lp = -out$value, which_opt = o)
  
  cat("finished optimization", o, "\n")
}

opt_samples = map(1:length(opts),
                  ~ ising_sweep(x, opts[[.]] %>% gather(name, b, 1:6) %>% pull(b) %>%
                                  setNames(names(b)), S, sample(10000000, 1)))

save.image(args$output)