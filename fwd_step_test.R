

source('fwd_step.R')
source('generate_data.R')
source('model_selection.R')

active_groups = function(groups, beta) {
  beta.ind = aggregate(beta, by=list(groups), FUN = function(beta.coords) any(beta.coords != 0))
  active.groups = beta.ind$Group.1[which(beta.ind$x == TRUE)]
  return(active.groups)
}


fwd_group_simulation = function(n, sigma, groups, beta, nsim, max.steps, alpha = .1, categorical = FALSE, rand.beta = FALSE, plot = FALSE) {

  weights = sqrt(rle(groups)$lengths)
  p = length(groups)
  g = length(unique(groups))
  true.active.groups = active_groups(groups, beta)
  num.nonzero = length(true.active.groups)
  upper = max(abs(beta))
  lower = min(abs(beta[beta != 0]))
  
  P.mat = matrix(nrow=nsim, ncol=max.steps)
  AS.mat = P.mat
  P.mat.b = P.mat
  AS.mat.b = P.mat
  recover.mat = P.mat
  mses.train = c()
  mses.test = c()
  mse.beta = c()
  if (rand.beta == FALSE) {
    # track fitted values
    beta.mat = matrix(0, nrow=nsim, ncol=p)
  }
  
  for (i in 1:nsim) {

    if (rand.beta == TRUE) {
      beta = beta_staircase(groups, num.nonzero, upper, lower,
        rand.sign = TRUE, permute = TRUE, perturb = TRUE)
      true.active.groups = active_groups(groups, beta)
    }

    if (categorical == TRUE) {
      X = categorical_design(n, groups, ortho.within = FALSE)
      X.test = categorical_design(n, groups, ortho.within = FALSE)
    } else {
      X = gaussian_design(n, groups)
      X.test = gaussian_design(n, groups)
    }

    Y = rnorm(n)*sigma
    Y.noiseless = X %*% beta
    Y.noiseless.test = X.test %*% beta
    Y.beta = Y.noiseless + Y
    Y.test = Y.noiseless.test + rnorm(n)*sigma

    # Null case
    results = forward_group(X, Y, groups, weights, sigma, max.steps)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    # Signal case
    results.b = forward_group(X, Y.beta, groups, weights, sigma, max.steps)
    P.mat.b[i, ] = results.b$p.vals
    AS.mat.b[i, ] = results.b$active.set
    recover.mat[i, ] = sapply(results.b$active.set, function(x)
                 is.element(x, true.active.groups))

    # Tracking prediction error and other things
    #
    #
    # this is broken, need to use stopping rules!
    #
    #
    groups.active = sapply(groups, function(x) is.element(x, results.b$active.set))
    fitted.model = lm(Y.beta ~ X[ , groups.active] - 1)
    fitted.beta = fitted.model$coefficients
    mses.train = c(mses.train, mean((fitted.model$fitted.values - Y.beta)^2))
    mses.test = c(mses.test, mean((predict(fitted.model, newdata = data.frame(X.test)) - Y.test)^2))
    mse.beta = c(mse.beta, sum((beta[groups.active] - fitted.beta)^2) + sum((beta[setdiff(1:p, which(groups.active))])^2))
    if (rand.beta == FALSE) {
      beta.mat[i, groups.active] = fitted.beta
    }
  }

  bar.quantiles <- c(.25, .75)
  point.quantiles <- c(.05, .95)
  null.Pvals <- colMeans(P.mat)
  null.Pvals.bar <- apply(P.mat, 2, function(col) quantile(col, probs = bar.quantiles))
  null.Pvals.point <- apply(P.mat, 2, function(col) quantile(col, probs = point.quantiles))

  Pvals = colMeans(P.mat.b)
  Pvals.bar <- apply(P.mat.b, 2, function(col) quantile(col, probs = bar.quantiles))
  Pvals.point <- apply(P.mat.b, 2, function(col) quantile(col, probs = point.quantiles))
  
  #recover.mat = apply(AS.mat.b, 1:2, function(x) is.element(x, true.active.groups))
  MSRS = colMeans(recover.mat)
  num.recovered.groups = rowSums(recover.mat[, 1:num.nonzero])
  MRR = mean(num.recovered.groups)
  print(MRR)

  if (plot == TRUE) {
    xax <- 1:max.steps
    nxax <- xax + 0.2
    plot.main <- paste("n, g, #nonzero, MRR =",
                       toString(paste(c(n, g, num.nonzero, MRR), sep = ", ")))
    if (rand.beta == TRUE) {
      plot.main = paste(plot.main, "(randomized signal)")
    }
    plot(xax, MSRS, type = "l", main = plot.main, xlab = "Step", ylab = "MSRS", ylim = c(-.1,1.1))
    abline(v = num.nonzero, lty = "dotted")
    
    points(nxax, null.Pvals, col="red")
    arrows(nxax, null.Pvals.bar[1, ], nxax, null.Pvals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "red")
    points(nxax, null.Pvals.point[1, ], col = "red", pch = 24, cex = .5)
    points(nxax, null.Pvals.point[2, ], col = "red", pch = 25, cex = .5)
  
    points(xax, Pvals, col="green")
    arrows(xax, Pvals.bar[1, ], xax, Pvals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "green")
    points(xax, Pvals.point[1, ], col = "green", pch = 24, cex = .5)
    points(xax, Pvals.point[2, ], col = "green", pch = 25, cex = .5)
  }

  if (rand.beta == FALSE) {
    return(list(null.p = P.mat, signal.p = P.mat.b, active.set = AS.mat.b, true.step = recover.mat, m1 = num.nonzero, mse.train = mses.train, mse.test = mses.test, mse.beta = mse.beta, beta.mat = beta.mat))
  } else {
    return(list(null.p = P.mat, signal.p = P.mat.b, active.set = AS.mat.b, true.step = recover.mat, m1 = num.nonzero, mse.train = mses.train, mse.test = mses.test, mse.beta = mse.beta))
  }
}

main = function() {

  nsim = 100

  # A small problem
  n = 50
  sigma = 1.01
  groups = c(1,1,2,2,3,4,4,4,5,5,6,7)
  upper = 0.9
  lower = 0.6
  num.nonzero = 3
  max.steps = 5
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE)
  garbage = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, plot = TRUE)

  # A larger problem
  n = 200
  sigma = 1.01
  groups = sort(c(rep(1:10, 2), rep(11:15, 5), rep(16:20, 10)))
  upper = 1
  lower = 0.2
  num.nonzero = 8
  max.steps = 10
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign = TRUE, permute = TRUE)
  garbage = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, plot = TRUE)

  # Randomized beta
  garbage = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
  
  # An n << p problem
  n = 100
  sigma = 1.01
  groups = sort(c(rep(1:100, 5), rep(101:150, 10)))
  upper = 3.4
  lower = 1.4
  num.nonzero = 10
  max.steps = 12
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign = TRUE)
  garbage = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, plot = TRUE)

}

if (interactive()) {
#  pdf('test_fwd_step.pdf')
#  main()
#  dev.off()
}

