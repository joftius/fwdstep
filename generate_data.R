#########################################
# Functions to generate simulation data #
#########################################

library(gtools)

SigmaSqrt = function(Sigma) {
  svdS = svd(Sigma)
  return(svdS$u %*% diag(sqrt(svdS$d)) %*% t(svdS$u))
}

# TODO
# * Random group sizes
true_active_groups = function(groups, beta) {
  beta.ind = aggregate(beta, by=list(groups), FUN = function(beta.coords) any(beta.coords != 0))
  active.groups = beta.ind$Group.1[which(beta.ind$x == TRUE)]
  return(active.groups)
}

# An approximation to start group lasso packages
lambda_max = function(X, Y, groups) {
  z = t(X) %*% Y
  return(max(aggregate(z, by=list(groups), FUN=function(x) sqrt(sum(x^2)) )$V1))
}

# Interfacing with grplasso package
grplasso_as = function(grplasso.fit, all.groups, tol=1e-10) {
  inds = which(abs(grplasso.fit$coefficients) > tol)
  return(unique(all.groups[inds]))
}

# active = true.active set, as = selected active set
ave_pow= function(num.nonzero, active, as) {
  return(length(intersect(active, as))/num.nonzero)
}

groupwise_center_scale = function(X, groups) {
  # Center by group
  ldata <- lapply(unique(groups), function(g) {
    inds <- which(groups == g)
    submat <- X[, inds]
    submat <- submat - mean(submat)
    return(submat)
  })
  X.out <- do.call(cbind, ldata)
  # Scale by group
  X.out = frob_normalize(X.out, groups)
  return(X.out)
}

# For coefficients of categorical variables,
# ensure zero-sum constraint.
center_rescale = function(v) {
    mod = sqrt(sum(v^2))
    if (mod > 0) {
        v = v - mean(v)
        newmod = sqrt(sum(v^2))
        if (newmod == 0) {
            stop("Coefficient is constant, cannot satisfy zero-sum")
        }
        v = v * mod / newmod
    }
    return(v)
}

constrain_beta = function(beta, groups, nz.groups, cat.groups) {
    group.labels = intersect(unique(nz.groups), cat.groups)
    for (g in group.labels) {
        gind = g == groups
        beta[gind] = center_rescale(beta[gind])
    }
    return(beta)
}

# Staircase signal
beta_staircase = function(groups, num.nonzero, upper, lower, rand.within=FALSE, rand.sign=TRUE, permute=TRUE, perturb=TRUE, cat.groups = NULL, staircase=TRUE) {
  # Generate a staircase-shaped signal vector
  # groups: vector of group indices in the form c(1,1,...,2,2,...)
  # num.nonzero: number of signal groups
  # upper: largest magnitude for a group
  # lower: lowest nonzero magnitude for a signal group
  # rand.within: if TRUE, adds small noise to coefficient for each group
  # rand.sign: if TRUE, randomly changes the signs of each coefficient
  # permute: if TRUE, permute indices of signal groups
  # perturb: if TRUE, adds small noise to all coefficients (making them nonconstant across group)
  labels <- unique(groups)
  g = length(labels)
  p = length(groups)

  beta = sort(seq(from=lower, to=upper, length=num.nonzero), decreasing=TRUE)

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

  nz.groups = labels[which(beta != 0)]
  # Expand
  beta = beta[groups]
  nz.inds = beta != 0

  if (rand.sign) {
    signs = sample(c(-1, 1), length(beta), replace = TRUE)
    beta = signs*beta
  }

  # Normalize coeff across group for fair comparison with non-grouped vars
  for (g in nz.groups) {
    group = g == groups
    gs = sum(group)
    beta[group] = beta[group]/sqrt(gs)
    bg.norm = sqrt(sum(beta[group]^2))

    # Add small noise in each group, but keep upper and lower fixed
    if (perturb) {
      beta[group] = beta[group] + rnorm(gs) * sqrt(1/10)
      bg.new.norm = sqrt(sum(beta[group]^2))
      beta[group] = beta[group] * bg.norm / bg.new.norm
    }
  }

  # Ensure coefs for categorical variables sum to 0
  # Maybe I don't need to do this
  if (length(cat.groups) > 0) {
      beta = constrain_beta(beta, groups, nz.groups, cat.groups)
  }

  if (!staircase) {
    #beta[nz.inds] = rnorm(sum(nz.inds))
    nnz = sum(nz.inds)
    beta[nz.inds] = sample(c(-1, 1), nnz, replace = TRUE)
    beta = beta/sqrt(nnz)
  }

  nz.labels <- labels[nz.groups]
  return(list(beta=beta, true.active=nz.labels, all.active=nz.labels))
}

# Normalize columns by 2-norm
col_normalize = function(X) {
  if (is.matrix(X)) {
    norms = sqrt(colSums(X^2))
    norms = ifelse(norms == 0, 1, norms)
    return(t(t(X) / norms))
  } else {
    norm = sqrt(sum(X^2))
    norm = ifelse(norm == 0, 1, norm)
    return(X/norm)
  }
}

# Normalize groups by Frobenius norm
frob_normalize = function(X, groups) {
    X.out = X
    for (g in unique(groups)) {
        inds = groups == g
        frob.norm = sqrt(sum(X[,inds]^2))
        if (is.na(frob.norm)) stop(paste("Matrix contains NA, check group", g))
        if (frob.norm > 0)
            X.out[,inds] = X.out[,inds]/frob.norm
    }
    return(X.out)
}

uniform_design = function(n, groups, col.normalize=FALSE, corr = 0) {
  p = length(groups)
  return(matrix(2*runif(n*p)-1, nrow=n))
}

# Fixed group sizes, gaussian design
gaussian_design = function(n, groups, col.normalize = FALSE, corr = 0) {
  p = length(groups)
  X = matrix(rnorm(n*p), nrow=n)

  ## if (ortho.within == TRUE) {
  ##   for (g in 1:max(groups)) {
  ##     group = groups == g
  ##     X[ , group] = svd(X[ , group])$u
  ##   }
  ## }

  if (corr != 0) {
    Z = matrix(rep(t(rnorm(n)), p), nrow=n)
    X = sqrt(1-corr)*X + sqrt(corr)*Z
  }

  if (col.normalize) {
    X = col_normalize(X)
  }
  return(X)
}

binary_design = function(n, groups) {
  p = length(groups)
  entries = sample(c(0,1), n*p, replace=TRUE)
  X = matrix(entries, nrow=n)
}

ternary_design = function(n, groups) {
  p = length(groups)
  entries = sample(c(-1,0,1), n*p, replace=TRUE)
  X = matrix(entries, nrow=n)
}

orthogonal_design = function(n, groups) {
  # If n > p, orthogonalize
  # If p = k*n, orthogonalize the k-submatrices
  # Else stop(error)
  X = gaussian_design(n, groups)
  if (n > p) {
    X = svd(X)$u
  } else {
    if (p %% n == 0) {
      k = p/n
      for (j in 0:(k-1)) {
        X[ , (j+1):(j+n)] = svd(X[ , (j+1):(j+n)])$u
      }
    } else {
      stop("p not a multiple of n")
    }
  }
  X = col_normalize(X)
  return(X)
}


cat_level_test = function(cat.levels) {

}
# Fixed group sizes, categorical design
# Important: binary requires two indices in groups, e.g. c(1,1,...)
categorical_design = function(n, groups, col.normalize = FALSE) {

  if (min(rle(groups)$lengths) <= 1) {
    stop("Minimum number of levels must be at least 2")
  }
  p = length(groups)
  X = matrix(nrow=n, ncol=p)
  X.cat = data.frame(ind=1:n)

  for (g in 1:max(groups)) {
    group = groups == g
    group.size = sum(group)
    cat.levels = rep(0, n)
    dir.prior = rep(0, group.size)
    # Sample from dirichlet prior
    dir.prior = rdirichlet(1, rep(group.size, group.size))
    # Resample until no levels are empty
    while ((min(table(cat.levels)) < 5) | (length(unique(cat.levels)) < group.size)) {
      cat.levels = sample(1:group.size, n,  replace=TRUE, prob = dir.prior)
    }
    X.cat = cbind(X.cat, factor(cat.levels))
    X[ , group] = unname(model.matrix(~ factor(cat.levels) - 1)[ , 1:group.size])
#    if (ortho.within == TRUE) {
#      X[ , group] = X[ , group] %*% diag(1/sqrt(colSums(X[ , group])))
#    }
  }
  X.cat = X.cat[,-1]
  names(X.cat) = paste0("X", 1:max(groups))
  if (col.normalize) {
    X = col_normalize(X)
  }
  #return(list(X=X, X.cat=X.cat))
  return(X)
}

# Input: inds of main effects
# Output: a permutation of those inds with no fixed points (derangement)
derangement = function(inds) {
  fixed = 1
  while (fixed > 0) {
    perm = sample(inds)
    fixed = sum(perm == inds)
  }
  return(perm)
}

