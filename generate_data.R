#########################################
# Functions to generate simulation data #
#########################################


# TODO
# * Random group sizes
# * Fix categorical data generation

# Staircase signal

beta_staircase = function(g, k, num.nonzero, upper, lower, rand.within=FALSE, rand.sign=FALSE) {
  # Generate a staircase-shaped signal vector
  # g: number of groups
  # k: group size
  # num.nonzero: number of signal groups
  # upper: largest magnitude for a group
  # lower: lowest nonzero magnitude for a signal group
  # rand.within: if TRUE, adds small noise to coefficients
  # rand.sign: if TRUE, randomly changes the signs of each coefficient
  q = k * num.nonzero
  p = k * g
  beta = sort(rep(seq(from=lower, to=upper, length=num.nonzero), k), decreasing=TRUE)

  if (rand.sign) {
    signs = sample(c(-1, 1), q, replace = TRUE)
    beta = signs*beta
  }

  if (rand.within) {
    mult = (upper - lower)/(num.nonzero - 1)
    mult = sqrt(mult / 6) # make noise small
    noise = rnorm(q) * mult
    beta = beta + noise
  }

  beta = c(beta, rep(0, p-q))
  return(beta)
}

# Fixed group sizes, gaussian design
simulate_fixed = function(n, g, k, orthonormal=TRUE, beta=0) {
  p = g*k
  X = matrix(rnorm(n*p), nrow=n)
  groups = list()
  for (group in 1:g) {
    groups[[group]] = k*(group-1) + 1:k
  }
  if (length(beta) == 1) {
    beta = matrix(0*1:p, ncol=1)
  }
  if (orthonormal == TRUE) {
    for (group in groups) {
      X[,group] = svd(X[,group])$u
    }
  }
  Y = X %*% beta + rnorm(n)
  weights = rep(1,g) * sqrt(k)
  return(list(X=X, Y=Y, groups=groups, weights=weights))
}


# Fixed group sizes, categorical design
simulate_fixed_cat = function(n, g, k, orthonormal=TRUE, beta=0) {
  p = g*k
  X = matrix(nrow=n, ncol=p)
  groups = list()
  for (group in 1:g) {
    groups[[group]] = k*(group-1) + 1:k
    cat.levels = 1
    # Resample until no levels are empty
    while (length(unique(cat.levels)) < k) {
      cat.levels = sample(1:k, n, replace=TRUE)
    }
    cat.binary = unname(model.matrix(~ factor(cat.levels) - 1)[1:n, 1:k])
    if (orthonormal == TRUE) {
      X[1:n, k*(group-1) + 1:k] = cat.binary %*% diag(1/sqrt(colSums(cat.binary)))
    } else {
      X[1:n, k*(group-1) + 1:k] = cat.binary
    }
  }
  if (length(beta) == 1) {
    beta = matrix(0*1:p, ncol=1)
  }
  Y = X %*% beta + rnorm(n)
  weights = rep(1,g) * sqrt(k)
  return(list(X=X, Y=Y, groups=groups, weights=weights))
}
