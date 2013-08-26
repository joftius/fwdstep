
trignometric_form = function(num, den, weight, tol=1.e-6) {
  
  a = num
  b = den
  w = weight
  
  norma = sqrt(sum(a^2))
  normb = sqrt(sum(b^2))
  if ((norma / normb) < tol) {
    return(c(0, Inf))
  }
  
  Ctheta = sum(a*b) / (norma*normb)
  Ctheta = min(max(Ctheta, -1), 1)
  Stheta = sqrt(1-Ctheta^2)
  theta = acos(Ctheta)
  
  Sphi = sqrt(sum(b^2)) * Stheta / w
          phi1 = asin(Sphi)
  phi2 = pi - phi1
  
  V1 = norma * cos(phi1) / (w - normb * cos(theta-phi1))
  V2 = norma * cos(phi2) / (w - normb * cos(theta-phi2))
  
  if (normb < w) {
    return(c(max(c(V1,V2)), Inf))
  }
  else {
    return(c(min(c(V1,V2)), max(c(V1,V2))))
  }
}

group_lasso_knot <- function(X, Y, groups, weights) {

  U = t(X) %*% Y
  g = length(weights)
  p = length(groups)
  terms = matrix(0, g)
                                        # first, compute the terms that go into
                                        # computing the dual norm
  for (i in 1:p) {
    terms[groups[i]] = terms[groups[i]] + U[i]^2
  }
  
  for (j in 1:g) {
    terms[j] = sqrt(terms[j]) / weights[j]
  }
  
  imax = which.max(terms)
  L = terms[imax]
  wmax = weights[imax]
  
  which = groups == imax
  
                                        # assuming rank == num of variables in group...
  kmax = sum(which)
  
  Uwhich = U[which]
  soln = (Uwhich / sqrt(sum(Uwhich^2))) / wmax
  
  Xmax = X[,which]
  
  if (kmax > 1) {
    soln = (Uwhich / sqrt(sum(Uwhich^2))) / wmax
    Xeta = (Xmax %*% soln)[,1]
    Xmax = Xmax - outer(Xeta, soln, '*') / sum(soln^2)
    Wmax = Xmax[,1:(ncol(Xmax)-1)]
    Xeta = lsfit(Wmax, Xeta, intercept=FALSE)$residuals
  }
  else {
    Xeta = Xmax / wmax * sign(U[which])
  }
  conditional_variance = sum(Xeta^2)
  Xeta = Xeta / conditional_variance
  
  C_X = t(X) %*% Xeta
  
  a = U - C_X * L
  b = C_X
  
  Vplus = c()
  Vminus = c()
  for (label in 1:g) {
    if (label != imax) {
      group = groups == label
      weight = weights[label]
      tf = trignometric_form(a[group], b[group], weight)
      Vplus = c(Vplus, tf[1])
      Vminus = c(Vminus, tf[2])
    }
  }
  
  if (length(Vplus) >= 1) {
    Vplus = Vplus[!is.nan(Vplus)]
    Vminus = Vminus[!is.nan(Vminus)]
    Mplus = max(Vplus)
    Mminus = min(Vminus)

  }
  else {
    Mplus = 0
    Mminus = Inf
  }
  return(list(L=L, Mplus=Mplus, Mminus=Mminus, var=conditional_variance, k=kmax))
}

pvalue <- function(L, Mplus, Mminus, sd, k, sigma=1) {
  first.term = pchisq((Mminus/(sd*sigma))^2, k, lower.tail=TRUE)
  if (first.term == 1) {
    num = pchisq((L/(sd*sigma))^2, k, lower.tail=FALSE, log.p=TRUE)
    den = pchisq((Mplus/(sd*sigma))^2, k, lower.tail=FALSE, log.p=TRUE)
    value = exp(num - den)
  } else {
    print("#### Something strange ####")
    #print(c(Mminus, first.term, L/(sd*sigma)^2, Mplus/(sd*sigma)^2))
    num = first.term - pchisq((L/(sd*sigma))^2, k, lower.tail=TRUE)
    den = first.term - pchisq((Mplus/(sd*sigma))^2, k, lower.tail=TRUE)
    value = num/den
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

