
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
generate_glint = function(X, groups, cat.groups = NULL, ...) {
  X.out = X
  n = dim(X)[1]
  g.out = groups
  main.groups = list()
  int.groups = list()
  main.out = rep(1, length(groups))
  gmax = max(groups)
  inds.out = matrix(rep(0, gmax^2), nrow=gmax)
  gnew = gmax
  for (g in 1:(gmax-1)) {
    ginds = which(groups == g)
    gs = length(ginds)
    for (h in (g+1):gmax) {
      gnew = gnew + 1
      int.groups[[g]] = gnew
      int.groups[[h]] = gnew
      main.groups[[gnew]] = c(g, h)
      Xgh = matrix(NA, nrow=n)
      hinds = which(groups == h)
      hs = length(hinds)
      gcont = ifelse(is.element(g, cat.groups), FALSE, TRUE)
      hcont = ifelse(is.element(h, cat.groups), FALSE, TRUE)

      if ((gcont) & (hcont)) {

        Xgh = cbind(Xgh, X[, c(ginds, hinds)])
        main.out = c(main.out, 0, rep(1, gs+hs), rep(0, gs*hs))
      } else if (gcont) {
        Xgh = cbind(Xgh, X[, hinds])
        main.out = c(main.out, rep(1, gs), rep(0, gs*hs))
      } else if (hcont) {
        Xgh = cbind(Xgh, X[, hinds])
        main.out = c(main.out, rep(1, hs), rep(0, gs*hs))
      } else {
        # Both are categorical
        main.out = c(main.out, rep(0, gs*hs))
      }
      for (i in ginds) {
        for (j in hinds) {
          Xij = X[, i] * X[, j]
          Xgh = cbind(Xgh, Xij)
        }
      }
      # Remove NA column used to initialize
      Xgh = Xgh[, -1]
      if ((gcont) & (hcont)) {
        Xgh = cbind(rep(1/n, n), Xgh)
      }
      X.out = cbind(X.out, Xgh)
      g.out = c(g.out, rep(gnew, ncol(Xgh)))
      inds.out[g, h] = gnew
      inds.out[h, g] = gnew
    }
  }
  colnames(X.out) = NULL
  return(list(X=X.out, main.groups=groups, all.groups=g.out, default.groups=main.groups, special.groups=inds.out, main.inds=main.out))
}

# Signal vector generation for glinternet
# convention: 1/3 of nonzero are only main effects, 2/3 have interactions
# k should be divisible by 3
# Note: at most (5/3)*k main effects and (2/3)*num.nz interactions
# are included, but group sparsity still = k
beta_glint = function(groups, all.groups, special.groups, k, num.default, upper, lower, rand.sign=TRUE, perturb=TRUE, cat.groups = NULL) {

  num.main = num.default
  int.groups = special.groups
  main.group.labels = unique(groups)
  main.inds = 1:length(groups)
  num.ints = k - num.main
  if (num.ints <= 0) {
    stop("At least 1 nonzero interaction is required")
  }

  p = dim(int.groups)[1]
  nz.mains = sample(main.group.labels, num.main, replace=FALSE)
  nz.main.inds = sapply(all.groups, function(x) is.element(x, nz.mains))
  remaining.mains = setdiff(main.group.labels, nz.mains)
  
  if (length(remaining.mains)/2 < num.ints) {
    stop("Not enough main effects to form desired number of interactions")
  }
  int.mains = sample(remaining.mains, 2*num.ints, replace=FALSE)

  nz.ints = c()
  nz.int.inds = rep(FALSE, length(all.groups))
  for (i in 1:num.ints) {
    this.int = glint_special_group_of(int.mains[i], int.mains[num.ints+i], int.groups)
    nz.ints = c(nz.ints, this.int)
    this.int.inds = which(this.int == all.groups)
    # Check if both main effects are continuous
    if (length(intersect(c(int.mains[i], int.mains[num.ints+i]), cat.groups)) == 0) {
      # In this case, don't add another coeff to the intercept
      this.int.inds = this.int.inds[-1]
    }
    nz.int.inds[this.int.inds] = TRUE
  }

  all.active = c(nz.mains, int.mains, nz.ints)
  nz.groups = sort(c(nz.mains, nz.ints))
  nz.inds = as.logical(nz.main.inds + nz.int.inds)
  if (length(nz.groups) != k) {
    stop("Incorrent number of nonzero groups, don't know why")
  }
  
  magnitudes = sample(seq(from=upper, to=lower, length=k))
  beta = rep(0, length(all.groups))

  for (i in 1:k) {
    ind.i = all.groups == nz.groups[i]
    beta[ind.i] = magnitudes[i]
  }
#  beta[nz.inds] = magnitudes
  
  if (rand.sign) {
    signs = sample(c(-1, 1), sum(nz.inds), replace = TRUE)
    beta[nz.inds] = signs*beta[nz.inds]
  }

  # Normalize coeff across group for fair comparison with non-grouped vars
  # Normalize interaction terms separately
  int.inds = main.inds != 1
  
  for (g in nz.groups) {
    group = g == all.groups
    if (g <= p) {
      gs = sum(group)
      beta[group] = beta[group]/sqrt(gs)
      bg.norm = sqrt(sum(beta[group]^2))
      
      if (perturb) {
        beta[group] = beta[group] + rnorm(gs) * sqrt(1/10)
        bg.new.norm = sqrt(sum(beta[group]^2))
        beta[group] = beta[group] * bg.norm / bg.new.norm
      }
      
    } else {
      group = nz.inds & group
      main.part = main.inds & group
      int.part = int.inds & group
      s.main = sum(main.part)
      s.int = sum(int.part)
      beta[main.part] = beta[main.part]/sqrt(s.main)
      # Inflate interaction terms to give them a chance
      beta[int.part] = sqrt(2)*beta[int.part]/sqrt(s.int)
      bg.main.norm = sqrt(sum(beta[main.part]^2))
      bg.int.norm = sqrt(sum(beta[int.part]^2))
      
      if (perturb) {
        beta[main.part] = beta[main.part] + rnorm(s.main) * sqrt(1/10)
        bg.main.new.norm = sqrt(sum(beta[main.part]^2))
        beta[main.part] = beta[main.part] * bg.main.norm / bg.main.new.norm
        
        beta[int.part] = beta[int.part] + rnorm(s.int) * sqrt(1/10)
        bg.int.new.norm = sqrt(sum(beta[int.part]^2))
        beta[int.part] = beta[int.part] * bg.int.norm / bg.int.new.norm
      }
      
    }
  }
  
##   # Normalize coeff across group for fair comparison with non-grouped vars
##   # Only do this if mu = X * beta is formed BEFORE normalizing X
##   for (g in unique(all.groups)) {
##     group = g == all.groups
##     gs = sum(group)
##     if (g <= p) {
##       # Main effects normalized as usual
##       beta[group] = beta[group]/sqrt(gs) #sqrt(max(1,gs-1))
##     } else {
##       # Interaction effects are normalized separately
##       # i.e. they do not share their coeff mass with the main effects
##       # otherwise it would be difficult to find them
##       mains.of.g = main_effects_of(g, int.groups)
##       mgs = 0
##       for (h in mains.of.g) {
##         mgs = mgs + sum(h == all.groups)
##       }
##       igs = gs - mgs
##       ginds = which(group)
##       start.ind = 1
##       # Check if both groups are continuous
##       if (length(intersect(mains.of.g, cat.groups)) == 0) {
##         # Don't add another coeff to the intercept
##         start.ind = 2
##       }
##       for (h in mains.of.g) {
##         hs = sum(h == all.groups)
##         end.ind = start.ind + hs - 1
##         beta[ginds[start.ind:end.ind]] = beta[ginds[start.ind:end.ind]]/sqrt(mgs)
##         start.ind = end.ind+1
##       }
##       # Inflate interaction size
##       beta[ginds[start.ind:gs]] = sqrt(2)*beta[ginds[start.ind:gs]]/sqrt(igs)
##     }
##   }

  return(list(beta=beta, true.ints=nz.ints, true.mains = nz.mains, true.active = nz.groups, all.active = all.active, true.special=int.mains))
}

#betad = beta_glinternet(all.groups, int.groups, k, upper, lower)
#beta = betad$beta
#true.ints = betad$true.ints
#for (x in true.ints) print(main_effects_of(x, int.groups))

# Find group index of interaction group for main effects h and g
glint_special_group_of = function(g, h, special.groups) {
  special.group = special.groups[g,h]
  return(special.group)
}

# Find group index of main effects for the interaction group gh
glint_default_group_of = function(gh, default.groups) {
  ## main.effs = which(special.groups == gh, arr.ind=TRUE)[1,]
  ## if (main.effs[1] < main.effs[2]) {
  ##   return(main.effs)
  ## } else {
  ##   return(rev(main.effs))
  ## }
  default.group = default.groups[[gh]]
  return(default.group)
}

