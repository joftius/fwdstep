
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
add_group = function(X, Y, groups, weights, Sigma, active.set = 0, eff.p = 0, cat.groups = NULL) {
  n = length(Y)
  # From group_lasso.R
  results = group_lasso_knot(X, Y, groups, weights, Sigma = Sigma, active.set = active.set)
  imax = results$i
  gmax = groups == imax
  kmax = results$k
  # Effective number of parameters
  new.eff.p = eff.p + sum(gmax)

  if (new.eff.p >= n) {
    stop("Too many variables added. Try decreasing max.steps")
  }

  # Also using group_lasso.R
  p.value = pvalue(results$L, results$Mplus, results$Mminus, sqrt(results$var), kmax)
  new.active.set = update_active_set(active.set, imax)

  # Form new residual
  X.project = X
  Xgmax = X[, gmax]
  Pgmax = Xgmax %*% ginv(Xgmax)
      
  # Project all other groups orthogonal to the one being added
  for (gind in 1:max(groups)) {
    if (gind != imax) {
      group = groups == gind
      X.project[, group] = X[, group] - Pgmax %*% X[, group]
    }
  }
  # If X[, gmax] is categorical, leave one column out  
  if ((is.element(imax, cat.groups)) & (is.matrix(Xgmax))) {
    Xgmax = Xgmax[ , -1]
    Pgmax = Xgmax %*% ginv(Xgmax)
  }
  Y.resid = lm(Y ~ Xgmax - 1)$residual
  #Y.resid = Y - Pgmax %*% Y
  
  return(list(test.output = results, var = results$var, p.value = p.value, added = imax, active.set = new.active.set, eff.p = new.eff.p, Y.update = Y.resid, X.update = X.project, grank=kmax))
}


# Iterate add_group for max.steps
forward_group = function(X, Y, groups, weights = 0, Sigma = NULL, max.steps = 0, cat.groups = NULL) {
  n = length(Y)
  group.sizes = rle(groups)$lengths

  if ((length(weights) == 1) & (weights[1] == 0)) {
    weights = sqrt(group.sizes)
  }

  # Estimate sigma instead?
  if (length(dim(Sigma)) == 0) {
    stop("Sigma needed here. Maybe try identity?")
  }
  
  if (max.steps == 0) {
    max.steps = length(unique(groups)) - 1
  }

  active.set = 0
  eff.p = 0
  p.vals = c()
  chi.pvals = c()
  c.vars = c()
  Ls = c()

  Y.update = Y
  X.update = X

  for (i in 1:max.steps) {
    output = add_group(X.update, Y.update, groups, weights, Sigma, active.set, eff.p, cat.groups = cat.groups)
    active.set = output$active.set
    imax = output$added
    grank = output$grank
    eff.p = output$eff.p
    RSS = sum(Y.update^2)
    Y.update = output$Y.update
    RSSdrop = RSS - sum(Y.update^2)
    chi.p = pchisq(RSSdrop, lower.tail=F, df=grank)
    #print(c(RSSdrop, round(grank), chi.p))
    X.update = output$X.update
    p.vals = c(p.vals, output$p.value)
    chi.pvals = c(chi.pvals, chi.p)
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

  return(list(active.set = active.set, p.vals = p.vals, chi.pvals = chi.pvals, Ls = Ls))
}

