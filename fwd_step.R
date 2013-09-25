
library(MASS) # for ginv
source('group_lasso.R')

update_active_set = function(active.set, group) {

  if (active.set[1] == 0) {
    new.active.set = group
  } else {
    new.active.set = c(active.set, group)
  }
  return(new.active.set)
  
}

add_group = function(X, Y, groups, weights, sigma, active.set = 0, eff.p = 0) {

  n = length(Y)
  results = group_lasso_knot(X, Y, groups, weights, active.set)
  imax = results$i
  gmax = groups == imax
  new.eff.p = eff.p + sum(gmax)

  if (new.eff.p >= n) {
    stop("Too many variables added. Try decreasing max.steps")
  }
  
  p.value = pvalue(results$L, results$Mplus, results$Mminus, sqrt(results$var), results$k, sigma=sigma)
  new.active.set = update_active_set(active.set, imax)

  Y.resid = lm(Y ~ X[, gmax])$residuals

  X.project = X
  gmax = groups == results$i
  Xgmax = X[, gmax]
  Pgmax = Xgmax %*% ginv(Xgmax)
  
  for (gind in 1:max(groups)) {
    if (gind != imax) {
      group = groups == gind
      X.project[, group] = (diag(rep(1, n)) - Pgmax) %*% X[, group] 
    }
  }

  return(list(test.output = results, p.value = p.value, added = imax, active.set = new.active.set, eff.p = new.eff.p, Y.update = Y.resid, X.update = X.project))
}



forward_group = function(X, Y, groups, weights = 0, sigma = 0, max.steps = 0) {

  n = length(Y)
  group.sizes = rle(groups)$lengths

  if ((length(weights) == 1) & (weights[1] == 0)) {
    weights = sqrt(rle(groups)$lengths)
  }

  # Estimate sigma instead?
  if (sigma == 0) {
    stop("Sigma estimate needed here")
  }
  
  if (max.steps == 0) {
    max.steps = length(unique(groups)) - 1
  }

  active.set = 0
  eff.p = 0
  p.vals = c()
  Ls = c()

  Y.update = Y
  X.update = X

  for (i in 1:max.steps) {
    
    output = add_group(X.update, Y.update, groups, weights, sigma, active.set, eff.p)
    active.set = output$active.set
    eff.p = output$eff.p
    Y.update = output$Y.update
    X.update = output$X.update
    p.vals = c(p.vals, output$p.value)
    # tracking lambda_2
    Ls = c(Ls, output$test.output[1])

    # Some overfitting considerations
    if ((eff.p >= n - max(group.sizes)) & (i < max.steps)) {
      if (eff.p > n - mean(group.sizes)) {
        warning("Overfitting imminent!")
      } else {
        warning("Overfitting may occur soon")
      }
    }
    
  }

  return(list(active.set = active.set, p.vals = p.vals, Ls = Ls))
}
