require(Rcpp)

sourceCpp("ising2.cpp")

ising_gibbs_derivs = function(x, mu, beta, S, seed) {
  NN = nrow(x) * ncol(x)
  res = ising_gibbs(x, mu, beta, S, seed)$states[(S / 2) : S,]
  
  f = matrix(0, nrow = ncol(res), ncol = 1)
  rownames(f) = colnames(res)
  jac = matrix(0, nrow = ncol(res), ncol = ncol(res) - 1)
  rownames(jac) = colnames(res)
  colnames(jac) = names(beta)
  for(i in 1:ncol(res)) {
    f[i, 1] = mean(res[, i]) / NN
    for(j in 2:ncol(res)) {
      jac[i, j - 1] = -cov(res[, i], res[, j]) / NN
    }
  }
  
  list(f = f, jac = jac)
}
