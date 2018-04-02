library(gmm)

setECIs(path, ecis)
runMC(path)
corrs = getCorrs(path)

Tfrac = getTemperatureFraction(path)

corrs %>%
  group_by(mu) %>%
  filter(mci > 500) %>%
  ungroup() %>%
  filter(mu == 8) %>%
  select(starts_with("corr")) %>%
  select(nonZero) %>%
  ggpairs

gfunc = function(b, x) {
  btmp = rep(0, length(ecis))
  btmp[nonZero] = b
  a = GgradG(btmp, getG = TRUE)
  cat("gfunc: ", b, "\n")
  t(a$g)
}

ggradGfunc = function(b, x) {
  btmp = rep(0, length(ecis))
  btmp[nonZero] = b
  a = GgradG(btmp, getG = TRUE)
  cat("ggradGfunc: ", b, "\n")
  a$Egrad[,nonZero]
}

gmm(gfunc,
    x = diag(rep(1, length(nonZero))),
    t0 = opts2[[1]][nonZero],
    gradv = ggradGfunc,
    type = "twoStep",
    method = "L-BFGS-B",
    traceIter = TRUE)

optim(opts2[[1]][nonZero], gfunc, ggradGfunc)

aa = runSimulation(ecis, getG = TRUE)
a = GgradG(ecis, getG = TRUE)

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

i = 1
setTemperatureFraction(path, ts[i])
aaa = GgradG2(ecis, data[[i]], getG = TRUE)
aa = GgradG(ecis, getG = TRUE)

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

getW = function(b) {
  g = GgradG(b, TRUE)$g
  w_inv = matrix(0, nrow = nrow(g), ncol = nrow(g))
  for(i in 1:ncol(g)) {
    w_inv = w_inv + g[, i] %*% t(g[, i])
  }
  w_inv = w_inv + 1e-8 * diag(nrow(g))
  solve(w_inv) / ncol(g)
  #w_inv
}
