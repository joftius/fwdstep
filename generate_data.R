#########################################
# Functions to generate simulation data #
#########################################


# TODO
# * Random group sizes
# * Fix categorical data generation
# *** beta coefs sum to 0 over categorical group!


# Staircase signal

beta_staircase = function(groups, num.nonzero, upper, lower, rand.within=FALSE, rand.sign=FALSE, permute=FALSE, perturb=FALSE, cat.vars = NULL) {
  # Generate a staircase-shaped signal vector
  # groups: vector of group indices in the form c(1,1,...,2,2,...)
  # num.nonzero: number of signal groups
  # upper: largest magnitude for a group
  # lower: lowest nonzero magnitude for a signal group
  # rand.within: if TRUE, adds small noise to coefficient for each group
  # rand.sign: if TRUE, randomly changes the signs of each coefficient
  # permute: if TRUE, permute indices of signal groups
  # perturb: if TRUE, adds small noise to all coefficients (making them nonconstant across group)
  g = max(groups)
  p = length(groups)
  beta = sort(seq(from=lower, to=upper, length=num.nonzero), decreasing=TRUE)

  if (rand.sign) {
    signs = sample(c(-1, 1), num.nonzero, replace = TRUE)
    beta = signs*beta
  }

  if (rand.within) {
    mult = (upper - lower)/(num.nonzero - 1)
    mult = sqrt(mult / 6) # make noise small
    noise = rnorm(num.nonzero) * mult
    beta = beta + noise
  }

  beta = c(beta, rep(0, g - num.nonzero))

  if (permute) {
    beta = sample(beta)
  }

  beta = beta[groups]
  nz.inds = beta != 0

  if (perturb) {
    beta[nz.inds] = beta[nz.inds] + rnorm(sum(nz.inds)) * sqrt(1/10)
    maxmod = max(abs(beta))
    beta[nz.inds] = beta[nz.inds] * upper / maxmod
  }

  # Ensure coefs for categorical variables sum to 0
  if (length(cat.vars) > 0) {

    for (g in cat.vars) {
      gind = groups == g & nz.inds
      if (sum(gind) > 0) {
        gmod = max(abs(beta[gind]))

        beta[gind] = beta[gind] - mean(beta[gind])
        gnewmod = max(abs(beta[gind]))

        if (gnewmod == 0) stop("Categorical variable with constant coeff (same for all levels)")
        beta[gind] = beta[gind] * gmod / gnewmod
      }
    }

  }

  return(beta)
}

# Fixed group sizes, gaussian design
gaussian_design = function(n, groups, ortho.within = FALSE) {
  p = length(groups)
  X = matrix(rnorm(n*p), nrow=n)

  if (ortho.within == TRUE) {
    for (g in 1:max(groups)) {
      group = groups == g
      X[ , group] = svd(X[ , group])$u
    }
  }

  X = X %*% diag(1/sqrt(colSums(X^2)))
  
  return(X)
}


# Fixed group sizes, categorical design
# Important: binary requires two indices in groups, e.g. c(1,1,...)
categorical_design = function(n, groups, ortho.within = FALSE) {

  if (min(rle(groups)$lengths) <= 1) {
    stop("Minimum number of levels must be at least 2")
  }
  p = length(groups)
  X = matrix(nrow=n, ncol=p)

  for (g in 1:max(groups)) {
    group = groups == g
    group.size = sum(group)
    cat.levels = NULL
    # Resample until no levels are empty
    while (length(unique(cat.levels)) < group.size) {
      cat.levels = sample(1:group.size, n, replace=TRUE)
    }
    X[ , group] = unname(model.matrix(~ factor(cat.levels) - 1)[ , 1:group.size])
#    if (ortho.within == TRUE) {
#      X[ , group] = X[ , group] %*% diag(1/sqrt(colSums(X[ , group])))
#    }

  }

  X = X %*% diag(1/sqrt(colSums(X^2)))
  
  return(X)
}

#n = 6
#groups = c(1,1,1,2,2,3,3,3,3)
#X = categorical_design(n, groups)

