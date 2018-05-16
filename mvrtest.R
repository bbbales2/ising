library(Matrix)
library(MASS)
library(tidyverse)

N = 10

Q = matrix(0, nrow = N, ncol = N)
for(i in 1:(N - 1)) {
  Q[i, i] = 30
  Q[i, i + 1] = -30
}
Q[N, N] = 30
Q[1, N] = 30
Q = Q + t(Q)
H = solve(Q * 1e0)

list(x = 1:N,
     y = mvrnorm(1, rep(0, N), H) > 0.0) %>%
  as.tibble %>%
  ggplot(aes(x, y)) +
  geom_point()

M = 10000
pc = rep(0, N)
for(m in 1:M) {
  pc = pc + (mvrnorm(1, rep(0, N), H) > 0.0) * 1
}
pc = pc / M

i2b = function(i) {
  intToBits(i)[1:N] %>%
    as.integer %>%
    as.vector
}

p = rep(0, N)
for(i in 1:(2^N)) {
  x = i2b(i)
  
  p[i] = exp(-t(x) %*% Q %*% x)
}
p = p / sum(p)

M = 10000
pd = rep(0, N)
for(m in 1:M) {
  pd = pd + i2b(sample(1:(2^N), 1, prob = p))
}
pd = pd / M

map(1:200, ~ i2b(sample(1:(2^N), 1, prob = p))) %>%
  do.call(rbind, .) %>%
  cov

map(1:200, ~ (mvrnorm(1, rep(0, N), H) > 0.0) * 1) %>%
  do.call(rbind, .) %>%
  cov

pd = rep(0, N)
for(i in 1:(2^N)) {
  x = intToBits(i)[1:N] %>%
    as.integer %>%
    as.vector
  
  pd = pd + p[i] * x
}

