require(Rcpp)

sourceCpp("ising2.cpp")

ising_gibbs_derivs = function(x, mu, beta, S, seed) {
  N = nrow(x) * ncol(x)
  ising_gibbs(x, mu, beta, S, seed)$states %>%
    as.tibble %>%
    summarize(dX0dQ1 = -cov(X0, Q1) / N,
              dX0dQ2 = -cov(X0, Q2) / N,
              dX0dQ3 = -cov(X0, Q3) / N,
              dX0dQ4 = -cov(X0, Q4) / N,
              dX0dQ5 = -cov(X0, Q5) / N,
              X0 = mean(X0) / N,
              Q1 = mean(Q1) / N,
              Q2 = mean(Q2) / N,
              Q3 = mean(Q3) / N,
              Q4 = mean(Q4) / N,
              Q5 = mean(Q5) / N)
}
