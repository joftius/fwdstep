
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
  cat.out = cat.groups
  gnew = gmax
#  omitted = 0
  for (g in 1:(gmax-1)) {
    ginds = which(groups == g)
    gs = length(ginds)
    for (h in (g+1):gmax) {
        main.append = c()

      Xgh = matrix(NA, nrow=n)
      hinds = which(groups == h)
      hs = length(hinds)
      gcont = !is.element(g, cat.groups)
      hcont = !is.element(h, cat.groups)
        gnew = gnew + 1
#      all.zero = FALSE

      if ((gcont) & (hcont)) {
        # Both continuous
        #Xgh = cbind(Xgh, rep(1,n))
        Xgh = cbind(Xgh, X[, c(ginds, hinds)])
        main.append = c(0, rep(1, gs+hs), rep(0, gs*hs))
      } else if ((gcont) & (!hcont)) {
        # g continuous
        Xgh = cbind(Xgh, X[, hinds])
        main.append = c(rep(1, hs), rep(0, gs*hs))
        cat.out = c(cat.out, gnew)
      } else if ((hcont) & (!gcont)) {
        # h continuous
        Xgh = cbind(Xgh, X[, ginds])
        main.append = c(rep(1, gs), rep(0, gs*hs))
        cat.out = c(cat.out, gnew)
      } else {
        # Both are categorical
        main.append = c(rep(0, gs*hs))
        cat.out = c(cat.out, gnew)
#        all.zero = TRUE
      }
        
      for (i in ginds) {
        for (j in hinds) {
          Xij = X[, i] * X[, j]
#          if (sum(Xij^2) != 0) {
#              all.zero = FALSE
#          }
          Xgh = cbind(Xgh, Xij)
        }
      }
      # Remove NA column used to initialize
      Xgh = Xgh[, -1]
#      if (all.zero) {
#        omitted = omitted + 1
#      } else {

          int.groups[[g]] = gnew
          int.groups[[h]] = gnew
          main.groups[[gnew]] = c(g, h)
          main.out = c(main.out, main.append)
# why is this here? Done above
#          if ((gcont) & (hcont)) {
#              Xgh = cbind(rep(1/n, n), Xgh)
#          }
          X.out = cbind(X.out, Xgh)
          g.out = c(g.out, rep(gnew, ncol(Xgh)))
          inds.out[g, h] = gnew
          inds.out[h, g] = gnew
#      }
    }
  }
  colnames(X.out) = NULL
  return(list(X=X.out, main.groups=groups, all.groups=g.out, default.groups=main.groups, special.groups=inds.out, main.inds=main.out)) #, omitted=omitted))
}


## for (i in 1:100) {
## b = beta_glint(groups, all.groups, special.groups, k, num.default, upper, lower, cat.groups = cat.groups)
## continc = is.element(11, b$all.active)
## print(paste("All cat?:", !continc))
## discreps = abs(sapply(unique(all.groups), function(g) sum(b$beta[g == all.groups])))
## failed = discreps > 1e-10
## print(paste("Constraints:", !any(failed)))
## if (any(failed)) print(which(failed))
## if (continc) print("-----------------------")
## }

# Signal vector generation for glinternet
# convention: 1/3 of nonzero are only main effects, 2/3 have interactions
# k should be divisible by 3
# Note: at most (5/3)*k main effects and (2/3)*num.nz interactions
# are included, but group sparsity still = k
beta_glint = function(groups, all.groups, special.groups, default.groups, k, num.default, upper, lower, rand.sign=TRUE, perturb=TRUE, cat.groups = NULL) {

  num.main = num.default
  int.groups = special.groups
  main.group.labels = unique(groups)
  main.inds = rep(FALSE, length(all.groups))
  main.inds[1:length(groups)] = TRUE
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
      g1 = int.mains[i]
      g2 = int.mains[num.ints+i]
      g1s = sum(groups == g1)
      g2s = sum(groups == g2)
    this.int = glint_special_group_of(g1, g2, int.groups)
    nz.ints = c(nz.ints, this.int)
    this.int.inds = which(this.int == all.groups)
    # Check if both main effects are continuous
    ints.cat = length(intersect(c(int.mains[i], int.mains[num.ints+i]), cat.groups))
    if (ints.cat == 0) {
        # In this case, don't add another coeff to the intercept
        #this.int.inds = this.int.inds[-1]
        # 
    } else if (ints.cat == 1) {
        # One of these is categorical
        #stop("Dunno yet how to constrain cont./cat. interaction")
    } else {
        # Both are categorical
        #this.int.inds = this.int.inds[(g1s+g2s+1):length(this.int.inds)]
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
  # Ensure categorical coeffs satisfy constraints
  int.inds = main.inds != 1
  for (g in nz.groups) {
      group = g == all.groups
      #print(c(length(group), length(nz.inds), length(main.inds)))
      if (g <= p) {
          gs = sum(group)
          beta[group] = beta[group]/sqrt(gs)
          bg.norm = sqrt(sum(beta[group]^2))
      
          if (perturb) {
              beta[group] = beta[group] + rnorm(gs) * sqrt(1/10)
              bg.new.norm = sqrt(sum(beta[group]^2))
              beta[group] = beta[group] * bg.norm / bg.new.norm
          }

          if (is.element(g, cat.groups)) {
              beta[group] = center_rescale(beta[group])
          }
      
      } else {
          #group = nz.inds & group
          main.part = main.inds & group
          int.part = int.inds & group
          s.main = sum(main.part)
          s.int = sum(int.part)
          if (s.main > 0) {
              beta[main.part] = beta[main.part]/sqrt(s.main)
              bg.main.norm = sqrt(sum(beta[main.part]^2))
          }
          # Inflate interaction terms to give them a chance
          beta[int.part] = sqrt(2)*beta[int.part]/sqrt(s.int)
          bg.int.norm = sqrt(sum(beta[int.part]^2))
      
          if (perturb) {
              if (s.main > 0) {
                  beta[main.part] = beta[main.part] + rnorm(s.main) * sqrt(1/10)
                  bg.main.new.norm = sqrt(sum(beta[main.part]^2))
                  beta[main.part] = beta[main.part] * bg.main.norm / bg.main.new.norm
              }
        
              beta[int.part] = beta[int.part] + rnorm(s.int) * sqrt(1/10)
              bg.int.new.norm = sqrt(sum(beta[int.part]^2))
              beta[int.part] = beta[int.part] * bg.int.norm / bg.int.new.norm
          }

          # Center categorical coefficients
          mains.of.g = glint_default_group_of(g, default.groups)
          cats.of.g = intersect(mains.of.g, cat.groups)
          if (length(cats.of.g) == 2) {
              # Remove global mean
              beta[int.part] = center_rescale(beta[int.part])
          } else if (length(cats.of.g) == 1) {
              # Which one to center? I think only the cat one,
              # i.e. first half of group
              cat.inds = which(group)
              cat.inds = cat.inds[1:(length(cat.inds)/2)]
              beta[cat.inds] = center_rescale(beta[cat.inds])
          }
          
      }  
  } # End for

  return(list(beta=beta, true.ints=nz.ints, true.mains = nz.mains, true.active = nz.groups, all.active = all.active, true.special=int.mains))
}

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

