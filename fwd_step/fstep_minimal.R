
# A minimal version of forward stepwise for phase diagram simulations

library(MASS)

# Append group to active set
update_active_set = function(active.set, group) {
  if (active.set[1] == 0) {
    new.active.set = group
  } else {
    new.active.set = c(active.set, group)
  }
  return(new.active.set)
}


# Determine which group to add next
add_group = function(X, R, groups, weights, active.set) {
  n = length(Y)
  g = length(weights)
  p = length(groups)
  
  U = t(X) %*% R
  terms = matrix(0, g)
  inner.prods = aggregate(U, by=list(groups), FUN=function(x) sum(x^2))$V1

  for (j in 1:g) {
    if (is.element(j, active.set)) {
      terms[j] = 0
    } else {
      terms[j] = inner.prods[j] / weights[j]
    }
  }
  imax = which.max(terms)
  return(imax)
}



# Input data, output sequence of added variables
fstep_fit = function(X, Y, groups, max.steps = 0) {
  n = length(Y)
  R = Y
  active.set = 0
  weights = sqrt(rle(groups)$lengths)
  
  if (max.steps == 0) {
    # Don't over-fit
    group.sizes = rle(groups)$values
    worst.case = sum(cumsum(sort(group.sizes, decreasing=TRUE)) <= n)
    max.steps = min(worst.case , length(unique(groups)))
  }

  for (step.count in 1:max.steps) {
    # Next group to add
    next.step = add_group(X, R, groups, weights, active.set)
    active.set = update_active_set(active.set, next.step)
    gmax = groups == next.step

    # Form new residual
    X.project = X
    Xgmax = X[, gmax]
    Xgmax.regress = Xgmax
    Pgmax = Xgmax %*% ginv(Xgmax)
    Hgmax = diag(rep(1, n)) - Pgmax
    # If X[, gmax] is categorical, leave one column out
    #if (sum(gmax) > 1) {
    #  if (length(unique(rowSums(Xgmax))) == 1) {
    #    Xgmax.regress = Xgmax[ , -1]
    #  }
    #}
    R = R - Pgmax %*% R

    ####### This is necessary for p-value? ########
    # Project all other groups orthogonal to the one being added
    for (gind in 1:max(groups)) {
      if (gind != next.step) {
        group = groups == gind
        X.project[, group] = Hgmax %*% X[, group]
      }
    }
    X.project = X.project %*% diag(1/sqrt(colSums(X.project^2)))
    X = X.project
  }

  return(active.set)
}


true_active_groups = function(groups, beta) {
  beta.ind = aggregate(beta, by=list(groups), FUN = function(beta.coords) any(beta.coords != 0))
  active.groups = beta.ind$Group.1[which(beta.ind$x == TRUE)]
  return(active.groups)
}
