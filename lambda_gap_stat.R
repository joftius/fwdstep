##############################################################
# Functions for computing test statistic from gap in lambdas #
##############################################################


# TODO
# * Check max_resid_proc, update for non-ortho / unequal weights
# * Check test_statistic
# * Is there numerical overflow? Why?


# Function to compute maximum of residual process wrt one group
max_resid_proc = function(Xh, wh, Xg, wg, Pg, y) {
  ug = t(Xg) %*% y
  ug = ug / sqrt(sum(ug^2))
  a = t(Xh) %*% (y - Pg %*% y) / wh
  Xug = Xg %*% ug
  b = (wg/wh) * (t(Xh) %*% Xug) / sum(Xug^2)
  ab = sum(a*b)
  norma2 = sum(a^2)
  normb2 = sum(b^2)
  # Rationalized denominator:
  return ((ab + sqrt(norma2*(1-normb2) + ab^2))/(1-normb2))
}


# Function to compute test statistic
test_statistic = function(X, Y, groups, weights, active.set=0) {
  inactive.groups = setdiff(1:length(groups), active.set)
  grad = t(X) %*% Y
  Ts = rep(-1,length(groups))
  for (i in inactive.groups) {
    Ts[i] = sqrt(sum(grad[groups[[i]]]^2)) / weights[i]
  }
  imax = which.max(Ts)
  gmax = groups[[imax]]
  X_gmax = X[,gmax]
  #P_gmax = X_gmax %*% ginv(X_gmax)     ### require(MASS)
  P_gmax = X_gmax %*% t(X_gmax)         ### Use orthogonality for now ###
  weight_gmax = weights[imax]
  Ms = c()
  for (i in inactive.groups) {
    if (i != imax) {
      Ms = c(Ms,max_resid_proc(X[,groups[[i]]], weights[i], X_gmax, weight_gmax, P_gmax, Y))
    }
  }
  M = max(Ms)
  u_gmax = grad[gmax]
  u_gmax = u_gmax / sqrt(sum(u_gmax^2))
  f_max = sum(u_gmax * grad[gmax]) / weight_gmax
  var_f_max = sum((X_gmax * u_gmax)^2) / weight_gmax^2
  T = f_max * (f_max - M) / var_f_max
  rank = sum(diag(P_gmax))
  return(list(T=T, lambda1=f_max, lambda2=M, gmax=gmax, imax=imax, weight=weight_gmax, rank=rank))
}


# Function to avoid numerical overflow
chisq_ratio = function(U, k) {
  numer.p = pchisq(U[1], k, lower.tail=FALSE, log.p=TRUE)
  denom.p = pchisq(U[2], k, lower.tail=FALSE, log.p=TRUE)
  return(exp(numer.p - denom.p))
}


# Function to compute p-value from test statistic
pval_from_lambdas = function(lambda1, lambda2, weight, rank) {
    U = weight*c(lambda1, lambda2)
    pval = chisq_ratio(U^2, rank)
  if (is.nan(pval)) {
    pval = 1
  }
  return(pval)
}

