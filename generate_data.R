#########################################
# Functions to generate simulation data #
#########################################


# TODO
# * Random group sizes
true_active_groups = function(groups, beta) {
  beta.ind = aggregate(beta, by=list(groups), FUN = function(beta.coords) any(beta.coords != 0))
  active.groups = beta.ind$Group.1[which(beta.ind$x == TRUE)]
  return(active.groups)
}

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

  if (rand.sign) {
    signs = sample(c(-1, 1), length(beta), replace = TRUE)
    beta = signs*beta
  }

  
  # Normalize coeff across group for fair comparison with non-grouped vars
  for (g in groups) {
    group = g == groups
    gs = sum(group)
    beta[group] = beta[group]/sqrt(gs)
  }

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

# Normalize columns by 2-norm
col_normalize = function(X) {
  norms = sqrt(colSums(X^2))
  return(t(t(X) / norms))
}

# Normalize groups by Frobenius norm
frob_normalize = function(X, groups) {
    X.out = X
    for (g in 1:max(groups)) {
        inds = groups == g
        frob.norm = sqrt(sum(X[,inds]^2))
        if (frob.norm > 0)
            X.out[,inds] = X.out[,inds]/frob.norm
    }
    return(X.out)
}

# Fixed group sizes, gaussian design
gaussian_design = function(n, groups, col.normalize = TRUE, corr = 0) {
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
  X = col_normalize(X)
}

ternary_design = function(n, groups) {
  p = length(groups)
  entries = sample(c(-1,0,1), n*p, replace=TRUE)
  X = matrix(entries, nrow=n)
  X = col_normalize(X)
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

  X = col_normalize(X)
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


#### Glinternet stuff
# X = [X1 X2 X3 ... X1X2 X1X3 X2X3 ...]
# all.groups = [1, 2, 3, ..., g+1, g+2, ...]
# main.groups = groups from original X (before interactions)
# int.groups = [[p+1, p+2,...], [...], ...]
#               [inds of interaction groups containing X1], [... X2], ...]
# beta nonzeros are for first k/3 main.groups


# Input design matrix with groups
# Output larger design matrix with all possible (grouped) interactions
# main.groups: groups of main effects (original variables)
# all.groups: include main effects and their copies with interactions
# int.groups: list, the [[g]]'th element contains all group indices that group g appears in
generate_glinternet = function(X, groups) {
  X.out = X
  n = dim(X)[1]
  g.out = groups
  main.out = groups
  gmax = max(groups)
  inds.out = matrix(rep(0, gmax^2), nrow=gmax)
  gnew = gmax
  for (g in 1:(gmax-1)) {
    for (h in (g+1):gmax) {
      ginds = which(groups == g)
      hinds = which(groups == h)
      gs = length(ginds)
      hs = length(hinds)
      X.out = cbind(X.out, X[, ginds]) #/sqrt(2))
      X.out = cbind(X.out, X[, hinds]) #/sqrt(2))
      gnew = gnew + 1
      for (i in ginds) {
        for (j in hinds) {
          Xij = X[, i] * X[, j]
          ## normij = sqrt(sum(Xij^2))
          ## if (normij > 0) {
          ##   Xij = Xij/normij
          ## }
          X.out = cbind(X.out, Xij)
        }
      }
      main.out = c(main.out, rep(g, gs), rep(h, hs), rep(0, gs*hs))
      g.out = c(g.out, rep(gnew, gs + hs + gs*hs))
      inds.out[g, h] = gnew
      inds.out[h, g] = gnew
    }
  }
  colnames(X.out) = NULL
  return(list(X=X.out, main.groups=groups, all.groups=g.out, int.groups=inds.out, main.inds=main.out))
}


# Signal vector generation for glinternet
# convention: 1/3 of nonzero are only main effects, 2/3 have interactions
# num.nonzero should be divisible by 3
# Note: at most (5/3)*num.nonzero main effects and (2/3)*num.nz interactions
# are included, but group sparsity still = num.nonzero
beta_glinternet = function(all.groups, int.groups, num.nonzero, upper, lower, rand.sign=TRUE, perturb=TRUE) {
  if (num.nonzero %% 3 != 0) stop("num.nonzero not divisible by 3")
  k = num.nonzero*(5/3)
  m = num.nonzero/3
  p = dim(int.groups)[1]
  magnitudes = seq(from=upper, to=lower, length=k)
  # Scramble the magnitudes around
  magnitudes = sample(magnitudes)
  true.ints = c()

  beta = rep(0, length(all.groups))
  beta[1:m] = magnitudes[1:m]
  for (i in (m+1):num.nonzero) {
    g = int_group_of(i, 2*m+i, int.groups)
    true.ints = c(true.ints, g)
    inds = all.groups == g
    beta[inds] = rep(magnitudes[i], sum(inds))
  }
  
  nz.inds = beta != 0

  if (rand.sign) {
    signs = sample(c(-1, 1), sum(nz.inds), replace = TRUE)
    beta[nz.inds] = signs*beta[nz.inds]
  }

  if (perturb) {
    beta[nz.inds] = beta[nz.inds] + rnorm(sum(nz.inds)) * sqrt(1/10)
    maxmod = max(abs(beta))
    beta[nz.inds] = beta[nz.inds] * upper / maxmod
  }
  
  # Normalize coeff across group for fair comparison with non-grouped vars
  # Don't do this, design matrix already used Frobenius norm
  for (g in unique(all.groups)) {
    group = g == all.groups
    gs = sum(group)
    if (g <= p) {
      # Main effects normalized as usual
      beta[group] = beta[group]/sqrt(gs) #sqrt(max(1,gs-1))
    } else {
      # Interaction effects are normalized separately
      # i.e. they do not share their coeff mass with the main effects
      # otherwise it would be difficult to find them
      mains.of.g = main_effects_of(g, int.groups)
      ginds = which(group)
      start.ind = 1
      for (h in mains.of.g) {
        hs = sum(h == all.groups)
        beta[ginds[start.ind:hs]] = beta[ginds[start.ind:hs]]/sqrt(hs)
        start.ind = start.ind + hs
      }
      # Inflate interaction size
      beta[ginds[start.ind:gs]] = sqrt(gs)*beta[ginds[start.ind:gs]] #/sqrt(gs-start.ind+1)
    }
  }

  return(list(beta=beta, true.ints=true.ints))
}

#betad = beta_glinternet(all.groups, int.groups, num.nonzero, upper, lower)
#beta = betad$beta
#true.ints = betad$true.ints
#for (x in true.ints) print(main_effects_of(x, int.groups))

# Find group index of interaction group for main effects h and g
int_group_of = function(g, h, int.groups) {
  int.group = int.groups[g,h]
  return(int.group)
}

# Find group index of main effects for the interaction group gh
main_effects_of = function(gh, int.groups) {
  main.effs = which(int.groups == gh, arr.ind=TRUE)[1,]
  if (main.effs[1] < main.effs[2]) {
    return(rev(main.effs))
  } else {
    return(main.effs)
  }
}

# input: g is group being added, p is # main effects
true_step_glinternet = function(g, p, int.groups, active.set, true.active.groups, already.counted) {
  true.active.int.groups = true.active.groups[which(true.active.groups > p)]
  if (g <= p) {
    # main effect group
    if (is.element(g, true.active.groups)) {
      # case 1
      if (!is.element(g, already.counted)) {
        return(1)
      }
    } else {
    # true active mixed groups that contain g
      containing.g = intersect(int.groups[g,], true.active.int.groups)
      if (length(containing.g) > 0) {
        if (!is.element(g, already.counted)) {
          # case 2
          return (1/3)
        }
      }
    }
    return(0)
  }
  mains.of.g = main_effects_of(g, int.groups)
  counted.mains = length(intersect(mains.of.g, already.counted))
  true.mains.of.g = length(intersect(mains.of.g, true.active.groups))
  if (is.element(g, true.active.groups)) {
    # case 4
    return(1-counted.mains/3)
  }
  if (true.mains.of.g > 0) {
    # case 5
    return(true.mains.of.g/3-counted.mains/3)
  }

  # other cases
  return(0)
}

# Truth: let m=k/3, first 1:m are main effect groups, next (m+1):2m main effects are matched to (2m+1):k with interactions
power_glinternet = function(k, true.ints, all.groups, int.groups, active.set) {
  m = k/3
  p = dim(int.groups)[1]
  S = 0
  counted.mains = c()
  for (g in active.set) {
    if (g <= k) {
      # Found a truly active main effect
      S = S + 1
      counted.mains = c(counted.mains, g)
    } else if (g > p) {
      
      if (is.element(g, true.ints)) {
        # Found a truly active interaction
        S = S + 1
      }
      mains.in.g = main_effects_of(g, int.groups)  
      for (h in mains.in.g) {
        if ((!is.element(h, counted.mains)) & (h <= k)) {
          # Also count the main effect for this interaction
          # if it hasn't already been counted
          S = S + 1
          counted.mains = c(counted.mains, h)
        }
      } 
    }
  }
  return(3*S/(4*k))
}
