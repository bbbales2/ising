data {
  int M; // Number of chemical potentials
  int P; // Number of parameters
  int N; // Number of samples in approx. distribution
  int Q; // Number of measured moments 
  real Tfrac;
  matrix[N, P] corrProb[M];
  matrix[Q, N] corrData[M];
  matrix[Q, P] corrDataEgrad[M];
  vector[Q] corrDataEg[M];
  vector[Q] cdata[M];
  vector[P] g0;
}

parameters {
  vector[P] g;
  real<lower = 0.0> sigma;
}

model {
  g ~ normal(0.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  
  for(m in 1:M) {
    //vector[N] p = softmax(-corrProb[m] * g / Tfrac);
    //cdata[m] ~ normal(corrData[m] * p + corrDataEgrad[m] * (g - g0), sigma);
    cdata[m] ~ normal(corrDataEg[m] + corrDataEgrad[m] * (g - g0), sigma);
    //rbind(t(p) %*% corrDataMats[[i]] - cdata[,i]))
  }
}