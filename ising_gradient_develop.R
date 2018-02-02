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

GgradG = function(g) {
  beta = g[1:5]
  gamma = g[6]

  a = mclapply(mus, function(mu) { smu(mu, beta, gamma) }, mc.cores = min(length(mus), 24))

  for(i in 1:length(a)) {
    J = length(a[[i]]$Eg)
    for(j in 1:J) {
      a[[i]]$g[j,] = a[[i]]$g[j,] - data[[i]]$f[j]
      a[[i]]$Eg[j] = a[[i]]$Eg[j] - data[[i]]$f[j]
    }
  }

  list(g = do.call("rbind", map(a, ~ .$g)),
       Eg = do.call("rbind", map(a, ~ .$Eg)),
       Egrad = do.call("rbind", map(a, ~ .$Egrad)))
}

getOmega = function(g) {
  omega = matrix(0, nrow = nrow(g$g), nrow(g$g))
  for(i in 1:ncol(g$g)) {
    omega = omega + g$g[,i] %*% t(g$g[,i])
  }
  omega / ncol(g$g) + 1e-6 * diag(nrow(omega))
}

#getOmega(g)

opts = list()
for(j in 11:40) {
  b = rnorm(6, 0.1, 0.25)
  eps = 0.1
  M = 500
  scaling = -log(0.01) / M
  for(i in 1:M) {
    g = GgradG(b)
    W = diag(nrow(g$g))
    u = t(g$Eg) %*% solve(W, g$Eg)
    grad = t(g$Eg) %*% solve(W, g$Egrad)
    dt = eps * exp(-i * scaling)
    
    cat("j : ", j, ", dt : ", dt, ", lp : ", u, ", params : ", b, "\n")
    b = b - dt * grad / (sqrt(sum(grad^2)) + 1e-10)
  }
  
  opts[[j]] = b
}
