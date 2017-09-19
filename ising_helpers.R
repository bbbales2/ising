require(Rcpp)

sourceCpp("ising.cpp")

ising_gibbs_derivs = function(x, mu, beta, gamma, kT, S, seed) {
  N = nrow(x)
  ising_gibbs(x, mu, beta, gamma, kT, S, seed)$states %>% as.tibble %>%
    summarize(mphi = mean(solo) / N^2,
              dmphi_dbeta = -cov(solo, pairs) / N^2,
              dmphi_dgamma = -cov(solo, corners) / N^2,
              d2mphi_dbeta2 = (mean(solo * pairs^2) +
                                 -mean(solo * pairs) * mean(pairs) +
                                 -cov(solo, pairs) * mean(pairs) +
                                 -var(pairs) * mean(solo)) / N^2,
              d2mphi_dgamma2 = (mean(solo * corners^2) +
                                  -mean(solo * corners) * mean(corners) +
                                  -cov(solo, corners) * mean(corners) +
                                  -var(corners) * mean(solo)) / N^2,
              d2mphi_dbetadgamma = (mean(solo * pairs * corners) +
                                      -mean(solo * pairs) * mean(corners) +
                                      -cov(solo, corners) * mean(pairs) +
                                      -cov(corners, pairs) * mean(solo)) / N^2)
}