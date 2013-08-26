
source('group_lasso.R')
source('generate_data.R')

pdf('test_group_lasso.pdf')
main = function() {

  n = 20
  p = 10
  sigma = 0.9
  groups = c(1,1,2,2,2,3,3,4,4,5)  
  weights = c(2,2.5,2,2,1.4)
  upper = 1
  lower = 0.999
  num.nonzero = 2
  beta = beta_staircase(4, 2, num.nonzero, upper, lower)
  beta = c(beta[1:3], beta[3:8], 0)
  print(beta)

  nsim = 1000
  P = c()
  P.beta = P
  for (i in 1:nsim) {
    X = matrix(rnorm(n*p),n,p)
    Y = rnorm(n)*sigma
    Y.beta = X %*% beta + Y

    results = group_lasso_knot(X, Y, groups, weights)
    P = c(P, pvalue(results$L, results$Mplus, results$Mminus, sqrt(results$var), results$k, sigma=sigma))

    results.b = group_lasso_knot(X, Y.beta, groups, weights)
    P.beta = c(P.beta, pvalue(results.b$L, results.b$Mplus, results.b$Mminus, sqrt(results.b$var), results.b$k, sigma=sigma))
  }

  qqplot(runif(nsim), P)
  abline(0,1)
  abline(h=0.1)
  qqplot(runif(nsim), P.beta)
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
    P = c(P, pvalue_MC(results$L, results$Mplus, results$Mminus,
      sqrt(results$var), results$k, sigma=sigma))
  }
  qqplot(runif(nsim), P)
  
}

#main_MC()

dev.off()
