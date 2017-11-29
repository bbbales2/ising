data {
  int N;
  int M;
  vector[N] X0;
  vector[M] Q;
  int useQ;
  real sigma;
}

parameters {
  vector[N] X0_;
  vector[M] Q_;
}

model {
  Q_ ~ normal(0, 1.0);
  X0_ ~ normal(X0, sigma);
  if(useQ > 0)
    Q_ ~ normal(Q, sigma);
}