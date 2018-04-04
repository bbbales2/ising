data {
  int M; // Number of chemical potentials
  int P; // Number of parameters
  int N; // Number of samples in approx. distribution
  int Q; // Number of measured moments 
  vector[P] corrDataEgrad[M, Q];
  real corrDataEg[M, Q];
  real cdata[M, Q];
  row_vector[P] g0;
}

parameters {
  row_vector[P] g;
  real<lower = 0.0> sigma;
}

model {
  g ~ normal(0.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  
  for(m in 1:M) {
    for(q in 1:Q) {
      cdata[m, q] ~ normal(corrDataEg[m, q] + (g - g0) * corrDataEgrad[m, q], sigma);// + 0.5 * to_row_vector(g - g0) * corrDataEgradgrad[m, q] * (g - g0)
    }
  }
}