########################################
# Functions for group-forward stepwise #
########################################

# TODO
# * Add function to do forward stepwise
# *** Proceed until a fixed stopping point (don't interpolate)
# *** Track order, pvals, gaps, etc, intelligently
# *** Decide how to orthogonalize & residualize
# * Check add1_group
# * Triple check everything


# Add the next group
add1_group = function(X, Y, groups, weights, active.set=0, residualize=TRUE, alpha=0.1, origX=X, origY=Y, betahat=0) {
  V = test_statistic(X, Y, groups, weights, active.set)
  V$X = X
  V$Y = Y
  V$betahat = betahat
  V$origX = origX
  V$origY = origY
  V$active.set = active.set
  imax = V$imax
  ginds = groups[[imax]]
  pval = pval_from_lambdas(V$lambda1, V$lambda2, V$weight, V$rank)
  V$pval = pval
  if ((length(betahat) > 1) & (length(active.set) < length(betahat) -1))  {
    V$active.set = setdiff(union(active.set, imax), 0)
    inds = unlist(groups[V$active.set])
    lmnew = lm(origY ~ origX[,inds] - 1)
    betanew = lmnew$coefficients
    V$betahat[inds] = betanew
    V$Y = lmnew$residuals
  }
  if (residualize == TRUE) {
    Xnew = X
    Xg = X[,ginds]
    Hg = diag(1,length(Y)) - Xg %*% t(Xg)
    for (g in setdiff(1:length(groups), imax)) {
      group = groups[[g]]
      Xnew[,group] = svd(Hg %*% X[,group])$u
    }
    V$X = Xnew
  }
  return(V)
}


