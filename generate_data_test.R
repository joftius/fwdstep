# Tests

source("generate_data.R")

main = function() {
  #n = 20
  #p = 10
  groups = c(1,1,2,2,2,3,3,4,4,5,6,6,6,6,7,8)
  #weights = c(2,2.5,2,2,1.4)
  num.nonzero = 4
  upper = 6
  lower = 5
  beta = beta_staircase(groups, num.nonzero, upper, lower)
  plot(beta, main = "Default")
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.within=TRUE)
  plot(beta, main = "Groups shifted")
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE)
  plot(beta, main = "Random sign")
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE, rand.within=TRUE)
  plot(beta, main = "Shifted and random sign")
  beta = beta_staircase(groups, num.nonzero, upper, lower, permute=TRUE)
  plot(beta, main = "Permuted")
  beta = beta_staircase(groups, num.nonzero, upper, lower, perturb=TRUE)
  plot(beta, main = "Perturbed")

  g = 8
  k = 5
  #data = simulate_fixed(100, g, k, beta=beta)
  #print(c(dim(data$X), length(data$Y), data$groups, data$weights))
}

pdf("test_generate_data.pdf")

main()

dev.off()
