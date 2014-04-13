
source('group_lasso.R')
source('generate_data.R')

pdf('test_group_lasso.pdf')
main = function() {

  n = 20
  sigma = 0.9
  groups = c(1,1,2,2,3,4,4,4,5,5,6,7)
  p = length(groups)
  weights = sqrt(rle(groups)$lengths)
  upper = 0.9
  lower = 0.4
  num.nonzero = 3
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE, permute=TRUE)
  print(beta)

  nsim = 1000
  P = c()
  P.beta = P
  signal.achieves = c()
  for (i in 1:nsim) {
    X = matrix(rnorm(n*p),n,p)
    Y = rnorm(n)*sigma
    Y.beta = X %*% beta + Y

    results = group_lasso_knot(X, Y, groups, weights)
    P = c(P, pvalue(results$L, results$lower_bound, results$upper_bound, sqrt(results$var), results$k, sigma=sigma))

    results.b = group_lasso_knot(X, Y.beta, groups, weights)
    group = groups == results.b$i
    signal.achieves = c(signal.achieves, any(beta[group] != 0))
    P.beta = c(P.beta, pvalue(results.b$L, results.b$lower_bound, results.b$upper_bound, sqrt(results.b$var), results.b$k, sigma=sigma))
  }

  P.small = P.beta <= 0.1
  outputs = c(mean(signal.achieves), mean(signal.achieves * P.small))
  qqplot(runif(nsim), P, main = "Null")
  abline(0,1)
  abline(h=0.1)

  qqplot(runif(nsim), P.beta, main = paste(c("Beta", unique(beta)), collapse = ", "),
         xlab = paste(c("Signal, detected", outputs), collapse =", "))
  abline(0,1)
  abline(h=0.1)
}

main()


main_MC = function() {

  n = 20
  p = 10
  sigma = 0.1

  nsim = 200
  P = c()
  for (i in 1:nsim) {
    X = matrix(rnorm(n*p),n,p)
    Y = rnorm(n)*sigma
    groups = c(1,1,2,2,2,3,3,4,4,5)
    weights = c(2,2.5,2,2,1.4)
    results = group_lasso_knot(X, Y, groups, weights)
    P = c(P, pvalue_MC(results$L, results$lower_bound, results$upper_bound,
      sqrt(results$var), results$k, sigma=sigma))
  }
  qqplot(runif(nsim), P)
  
}

#main_MC()

dev.off()
