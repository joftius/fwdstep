# Does this need to be modified for glinternet?

# This file contains core functions for computing the test statistic

library(Matrix)
library(MASS)

trignometric_form = function(num, den, weight, tol=1.e-10) {
  
  a = num
  b = den
  w = weight
  norma = sqrt(sum(a^2))
  normb = sqrt(sum(b^2))

  if (normb == 0) {
#    stop("Something is wrong, norm(b) shouldn't be zero!")
    return(c(0, Inf))
  }
  
  if ((norma / normb) < tol) {
    return(c(0, Inf))
  }
  
  Ctheta = sum(a*b) / (norma*normb)
  # Why was this next line here?
  #Ctheta = min(max(Ctheta, -1), 1)
  Stheta = sqrt(1-Ctheta^2)
  theta = acos(Ctheta)
  
  Sphi = normb * Stheta / w

  if (Sphi > 1) {
    # infeasible
    return(c(0, Inf))
  }
  
  phi1 = asin(Sphi)
  phi2 = pi - phi1
  
  V1 = norma * cos(phi1) / (w - normb * cos(theta-phi1))
  V2 = norma * cos(phi2) / (w - normb * cos(theta-phi2))
  
  if (normb < w) {
    # encode this infeasibility as Inf
    return(c(max(c(V1,V2)), Inf))
  }
  else {
    return(c(min(c(V1,V2)), max(c(V1,V2))))
  }
}

group_lasso_knot <- function(X, Y, groups, weights, Sigma, active.set=0) {

  U = t(X) %*% Y
  g.labels = unique(groups)
  g = length(weights)
  p = length(groups)
  # first, compute the terms that go into
  # computing the dual norm
  terms = sapply(g.labels, function(x) sqrt(sum(U[groups==x]^2))/weights[x])
  active.terms = sapply(g.labels, function(x) ifelse(is.element(x, active.set), 0, 1))
  terms = terms*active.terms
##   terms = matrix(0, g)  
##   for (i in 1:p) {
##     terms[groups[i]] = terms[groups[i]] + U[i]^2
##   }
##   for (j in 1:g) {
##     if (is.element(j, active.set)) {
##       terms[j] = 0
##     } else {
##       terms[j] = sqrt(terms[j]) / weights[j]
##     }
##   }
  imax = which.max(terms)
  L = terms[imax]
  if (L <= 0) {
    stop("Lambda should not be zero")
  }
    
  wmax = weights[imax]
  which = groups == imax
  Uwhich = U[which]
  Xmax = X[,which]
  kmax = sum(which)
  kmaxrank = kmax
  if (length(dim(Xmax)) == 2) {
    kmaxrank = rankMatrix(Xmax)[1]
  }

  eta = rep(0, p)
  eta[which] = (Uwhich / sqrt(sum(Uwhich^2))) / wmax
  if (kmax > 1) {
    # P = ....
    tangent.space = Null(Uwhich)
    # dim(tangent.space) = kmax \times (kmax-1)
    if (ncol(tangent.space) != kmax-1) {
      stop(paste0("Tanget space dimension is wrong: ", kmax,
                  ncol(tangent.space)))
    }
    tangent.space = cbind(tangent.space, rep(0, kmax))
    # Need ncol(V) = p
    V = matrix(0, ncol=kmax, nrow=p)
    V[which,] = tangent.space
    XV = X %*% V[,-ncol(V)]
    XV = cbind(XV, rep(0, nrow(XV)))
    XXy = Xmax %*% Uwhich
    if (!is.matrix(Sigma)) {
      SXV = Sigma*XV
      P = SXV %*% ginv(t(XV) %*% SXV) %*% t(XV)
      OP = diag(rep(1, ncol(P))) - P
      Hg = OP*Sigma
    } else {
      SXV = Sigma %*% XV
      P = SXV %*% ginv(t(XV) %*% SXV) %*% t(XV)
      OP = diag(rep(1, ncol(P))) - P
      Hg = OP %*% Sigma
    }
    Xeta = Hg %*% XXy / (wmax * sqrt(sum(Uwhich^2)))
    conditional_variance = t(XXy) %*% Hg %*% XXy
    # diag() coerces 1x1 matrix to numeric
    conditional_variance = diag(conditional_variance) / (wmax^2*sum(Uwhich^2))
  } else {
    # P = 0
    Xeta = Xmax / wmax * sign(Uwhich)
    conditional_variance = sum(Xeta^2)
  }
  # use formula above display (42) in tests:adaptive
  Xeta = Xeta / conditional_variance
  C_X = t(X) %*% Xeta

  a = U - C_X * L
  b = C_X
  
  Vplus = c()
  Vminus = c()
  nm.a = c()
  nm.b = c()
  nm.labels = setdiff(c(1:g), c(active.set, imax))
  for (label in nm.labels) {
    group = groups == label
    weight = weights[label]
    nm.a = c(nm.a, sqrt(sum(a[group]^2)))
    nm.b = c(nm.b, sqrt(sum(b[group]^2)))
    tf = trignometric_form(a[group], b[group], weight)
    Vplus = c(Vplus, tf[1])
    Vminus = c(Vminus, tf[2])
  }
#  print(rbind(nm.labels, nm.a, nm.b))
  
  if (length(Vplus) >= 1) {
    Vplus = Vplus[!is.nan(Vplus)]
    Vminus = Vminus[!is.nan(Vminus)]
    Mplus = max(Vplus)
    Mminus = min(Vminus)
  } else {
    Mplus = 0
    Mminus = Inf
  }
  return(list(L=L, Mplus=Mplus, Mminus=Mminus, var=conditional_variance, k=kmaxrank, i=imax))
}

# why using (sd*sigma) ???
pvalue <- function(L, Mplus, Mminus, sd, k) {
  first.term = pchisq((Mminus/sd)^2, k, lower.tail=TRUE)
  if (first.term == 1) {
    num = pchisq((L/sd)^2, k, lower.tail=FALSE, log.p=TRUE)
    den = pchisq((Mplus/sd)^2, k, lower.tail=FALSE, log.p=TRUE)
    value = exp(num - den)
  } else {
    #print(c(Mminus, first.term, L/(sd*sigma)^2, Mplus/(sd*sigma)^2))
    num = first.term - pchisq((L/sd)^2, k, lower.tail=TRUE)
    den = first.term - pchisq((Mplus/sd)^2, k, lower.tail=TRUE)
    value = num/den
    #print(paste("Mminus, first.term, value:", Mminus, first.term, value))
  }
  return(value)
}



q_0 <- function(M, Mminus, k, nsim=100) {
  Z = abs(rnorm(nsim))
  keep = Z < (Mminus - M)
  proportion = mean(keep)
  Z = Z[keep]
  exponent = - M*Z - M^2/2.
  if (k > 1) {
    exponent = exponent + (k-1)*log(Z+M)
  }
  
  C = max(exponent)
  
  return(list(E=mean(exp(exponent - C)) * proportion, C=C))
}

Q_0 = function(L, Mplus, Mminus, H, nsim=100) {
  
  result1 = q_0(L, Mminus, H, nsim=nsim)
  result2 = q_0(Mplus, Mminus, H, nsim=nsim)
  
  return(exp(result1$C-result2$C) * result1$E / result2$E)
}


pvalue_MC <- function(L, Mplus, Mminus, sd, k, sigma=1, nsim=50000) {
  sd = sd * sigma
  if (k > 1) {
    return(Q_0(L / sd, Mplus / sd, Mminus / sd, k, nsim=nsim))
  }
  else {
    return(Q_0(L / sd, Mplus / sd, Mminus / sd, 1, nsim=nsim))
  }
}

