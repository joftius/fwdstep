
# Solve x/log(x) = C
newton_solver = function(x0, C, max.steps=5) {
  x = x0
  for (i in 1:max.steps) {
    x = x - (x-C*log(x))/(1-1/log(x))
  }
  return(x)
}

mu_to_n = function(mu) {
  x0 = (4/mu)^2
  C = (2/mu)^2
  n = newton_solver(x0, C)
  return(round(n))
}

n_to_mu = function(n) {
  return(2*sqrt(log(n)/n))
}

mu_to_k = function(mu) {
  return(floor((1/mu+1)/2))
}

mu.max = .43
mu.min = .23
wanted.mus = seq(from = mu.max, to = mu.min, length = 14)

ns = mu_to_n(wanted.mus)
mus = n_to_mu(ns)
ks = mu_to_k(mus)
