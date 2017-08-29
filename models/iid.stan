data {
  int N;
  vector[N] y;
  real sigma;
}

parameters {
  vector[N] yhat;
}

model {
  y ~ normal(yhat, sigma);
}