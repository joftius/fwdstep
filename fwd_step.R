
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
  new.active.set = update_active_set(active.set, imax)
  inactive.set = setdiff(unique(groups), new.active.set)
  # Effective number of parameters
  new.eff.p = eff.p + sum(gmax)

  if (new.eff.p >= n) {
    stop("Too many variables added. Try decreasing max.steps")
  }

  conditional_variance = results$var

  # Form new residual
  X.project = X
  Xgmax = X[, gmax]
  Pgmax = Xgmax %*% ginv(Xgmax)
      
  # Project all other groups orthogonal to the one being added
  for (gind in inactive.set) {
    group = groups == gind
    X.project[, group] = X[, group] - Pgmax %*% X[, group]
  }

  Y.resid = Y - Pgmax %*% Y

  if (conditional_variance <= 0) {
    p.value = 1
  } else {
    p.value = pvalue(results$L, results$lower_bound, results$upper_bound, sqrt(conditional_variance), kmax)
  }
  
  return(list(test.output = results, var = results$var, p.value = p.value, added = imax, active.set = new.active.set, eff.p = new.eff.p, Y.update = Y.resid, X.update = X.project, grank=kmax))
}


# Bare bones version for power simulation without p-values
add_group_nopval= function(X, Y, groups, weights, Sigma, active.set = 0, eff.p = 0, cat.groups = NULL) {
    n = length(Y)
    g.labels = unique(groups)
    inactive.groups = setdiff(g.labels, active.set)

    terms = sapply(inactive.groups, function(x) sqrt(sum((t(X[,groups==x]) %*% Y)^2))/weights[x])

      imax = inactive.groups[which.max(terms)]
    
    gmax = groups == imax
    grank = sum(gmax)

    new.active.set = update_active_set(active.set, imax)
    inactive.set = setdiff(inactive.groups, imax)
    
    # Effective number of parameters
    new.eff.p = eff.p + sum(gmax)

    if (new.eff.p >= n) {
        stop("Too many variables added. Try decreasing max.steps")
    }

    # Form new residual
    X.project = X
    Xgmax = X[, gmax]
    Pgmax = Xgmax %*% ginv(Xgmax)
      
    # Project all other groups orthogonal to the one being added
    for (gind in inactive.set) {
        group = groups == gind
        X.project[, group] = X[, group] - Pgmax %*% X[, group]
    }

    Y.resid = Y - Pgmax %*% Y
  
    return(list(test.output = NA, var = NA, p.value = NA, added = imax, active.set = new.active.set, eff.p = new.eff.p, Y.update = Y.resid, X.update = X.project, grank=grank))
}


# Main forward stepwise function
# Iterate add_group for max.steps
forward_group = function(X, Y, groups, weights = 0, Sigma = NULL, max.steps = 0, cat.groups = NULL, tol=1e-10, pval=TRUE) {
  n = length(Y)
  group.sizes = rle(groups)$lengths

  if (length(weights) == 1) {
      print(weights)
      if (weights == 0) {
          weights = sqrt(group.sizes)
      } else if (weights == 1)  {
          weights = rep(1, length(group.sizes))
      } else {
          stop("Misspecified weights")
      }
  }

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
  ranks = c()
  Ynorms = c()

  Y.update = Y
  X.update = X

  for (i in 1:max.steps) {
      if (pval) {
          output = add_group(X.update, Y.update, groups, weights, Sigma, active.set, eff.p, cat.groups = cat.groups)

      } else {
          output = add_group_nopval(X.update, Y.update, groups, weights, Sigma, active.set, eff.p, cat.groups = cat.groups)
      }
      
      active.set = output$active.set
      imax = output$added
      grank = output$grank
      ranks = c(ranks, grank)
      eff.p = output$eff.p
      RSS = sum(Y.update^2)
      Y.update = output$Y.update
      X.update = output$X.update
      Ynorms = c(Ynorms, RSS)    
      RSSdrop = RSS - sum(Y.update^2)
      
      if (pval) {
          p.vals = c(p.vals, output$p.value)
          chi.p = pchisq(RSSdrop, lower.tail=F, df=grank)
          chi.pvals = c(chi.pvals, chi.p)
          c.vars = c(c.vars, output$var)
          Ls = c(Ls, output$test.output[[1]])
      }


      # Early stopping
      if ((RSSdrop < tol) & (i < max.steps)) {
          left = max.steps - i
          active.set = c(active.set, rep(0, left))
          p.vals = c(p.vals, rep(NA, left))
          chi.pvals = c(chi.pvals, rep(NA, left))
          Ls = c(Ls, rep(0, left))
          return(list(active.set = active.set, p.vals = p.vals, chi.pvals = chi.pvals, Ls = Ls, RSS = Ynorms))
      }

      # Some overfitting considerations
      if ((eff.p >= n - max(group.sizes)) & (i < max.steps)) {
          if (eff.p > n - mean(group.sizes)) {
              warning("Overfitting imminent!")
          } else {
              warning("Overfitting may occur soon")
          }
      }
  }

  outlist = list(active.set = active.set, p.vals = p.vals, chi.pvals = chi.pvals, Ls = Ls, RSS = Ynorms)

  return(outlist)
}

step.as = function(X, Y, steps, k) {
  df = data.frame(X)
  full.model = lm(Y ~ .-1, data=df)
  step.fit = step(lm(Y ~ 0, data=df),  scope=list(lower=~0, upper=full.model), direction="forward", k = k, steps = steps, trace = 0)
  step.as = as.numeric(gsub("X", "", names(step.fit$coefficients), fixed=TRUE))
  return(step.as)
}
