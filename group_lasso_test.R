
require('group_lasso.R')


pdf('group_lasso.pdf')
main = function() {

  n = 20
  p = 10
  sigma = 0.1
  
  nsim = 10000
  P = c()
  for (i in 1:nsim) {
    X = matrix(rnorm(n*p),n,p)
    Y = rnorm(n)*sigma
    groups = c(1,1,2,2,2,3,3,4,4,5)
    weights = c(2,2.5,2,2,1.4)
    results = group_lasso_knot(X, Y, groups, weights)
    P = c(P, pvalue(results$L, results$Mplus, results$Mminus, sqrt(results$var), results$k, sigma=sigma))
  }
  qqplot(runif(nsim), P)
  
}

main()


main_MC = function() {

  n = 20
  p = 10
  sigma = 0.1

  nsim = 1000
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

main_MC()


dev.off()
