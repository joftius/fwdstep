
library(MASS) # for ginv
source('group_lasso.R')

# Add a group to the active set
update_active_set = function(active.set, group) {
  if (active.set[1] == 0) {
    new.active.set = group
  } else {
    new.active.set = c(active.set, group)
  }
  return(new.active.set)
}

# Indicator of which groups are in the active set
active_groups = function(active.set, groups) {
  return(sapply(groups, function(x) is.element(x, active.set)))
}

# Indicator (or list?) of groups with at least one nonzero coefficient in beta
true_active_groups = function(groups, beta) {
  beta.ind = aggregate(beta, by=list(groups), FUN = function(beta.coords) any(beta.coords != 0))
  active.groups = beta.ind$Group.1[which(beta.ind$x == TRUE)]
  return(active.groups)
}

# Compute which group to add to the active set and associated pvalue
add_group = function(X.orig, X, Y, groups, weights, sigma, active.set = 0, eff.p = 0) {
  n = length(Y)
  # From group_lasso.R
  results = group_lasso_knot(X, Y, groups, weights, active.set = active.set)
  imax = results$i
  gmax = groups == imax
  # Effective number of parameters
  new.eff.p = eff.p + sum(gmax)

  if (new.eff.p >= n) {
    stop("Too many variables added. Try decreasing max.steps")
  }

  # Also using group_lasso.R
  p.value = pvalue(results$L, results$Mplus, results$Mminus, sqrt(results$var), results$k, sigma=sigma)
  new.active.set = update_active_set(active.set, imax)

  # Form new residual
  X.project = X
  Xgmax = X[, gmax]
  Xgmax.regress = Xgmax
  # If X[, gmax] is categorical, leave one column out
  ########### Is this doing anything? Xgmax.regress isn't used ###########
  if (sum(gmax) > 1) {
    if (length(unique(rowSums(Xgmax == 0))) == 1) {
#      print("categorical variable")
      Xgmax.regress = Xgmax[ , -1]
    }
  }
  Pgmax = Xgmax.regress %*% ginv(Xgmax.regress)
  Y.resid = (diag(rep(1,n)) - Pgmax) %*% Y
      
####### This is necessary for p-value? ########
  # Project all other groups orthogonal to the one being added
  for (gind in 1:max(groups)) {
    if (gind != imax) {
      group = groups == gind
      X.project[, group] = (diag(rep(1, n)) - Pgmax) %*% X[, group]
    }
  }
  # Renormalize
#  X.project = X.project %*% diag(1/sqrt(colSums(X.project^2)))
  
  return(list(test.output = results, var = results$var, p.value = p.value, added = imax, active.set = new.active.set, eff.p = new.eff.p, Y.update = Y.resid, X.update = X.project))
}


# Iterate add_group for max.steps
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
  c.vars = c()
  Ls = c()

  Y.update = Y
  X.update = X
  X.orig = X

  for (i in 1:max.steps) {
    
    output = add_group(X.orig, X.update, Y.update, groups, weights, sigma, active.set, eff.p)
    active.set = output$active.set
    eff.p = output$eff.p
    Y.update = output$Y.update
    X.update = output$X.update
    p.vals = c(p.vals, output$p.value)
    c.vars = c(c.vars, output$var)
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

