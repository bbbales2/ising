data {
  int N;
  real x[N];
  vector[N] y;
  real sigma;
  real l;
}

parameters {
  vector[N] yhat;
}

model {
  y ~ multi_normal(yhat, cov_exp_quad(x, sigma, l));
}