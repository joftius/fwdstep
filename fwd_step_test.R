
source('fwd_step.R')
source('generate_data.R')

pdf('test_fwd_step.pdf')
main = function() {

  n = 20
  sigma = 0.9
  groups = c(1,1,2,2,3,4,4,4,5,5,6,7)
  p = length(groups)
  max.steps = 5
  weights = sqrt(rle(groups)$lengths)
  upper = 0.9
  lower = 0.4
  num.nonzero = 3
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE)
  print(beta)

  nsim = 200
  P.mat = matrix(nsim, max.steps)
  AS.mat = P.mat
  P.mat.b = P.mat
  AS.mat.b = P.mat
  
  for (i in 1:nsim) {
    X = matrix(rnorm(n*p),n,p)
    Y = rnorm(n)*sigma
    Y.beta = X %*% beta + Y

    results = forward_group(X, Y, groups, weights, max.steps)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    results.b = forward_group(X, Y.beta, groups, weights, max.steps)
    P.mat.b[i, ] = results.b$p.vals
    AS.mat.b[i, ] = results.b$active.set

  }
  
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
