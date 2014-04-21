#########################################
# Functions to generate simulation data #
#########################################

library(splines)
source('generate_data.R')


#### GAMsel stuff
# X = [X1 X2 X3 ... X1X2 X1X3 X2X3 ...]
# all.groups = [1, 2, 3, ..., g+1, g+2, ...]
# main.groups = groups from original X (before interactions)
# int.groups = [[p+1, p+2,...], [...], ...]
#               [inds of interaction groups containing X1], [... X2], ...]
# beta nonzeros are for first k/3 main.groups


##### This only supports groups of size 1 right now #####
# Input design matrix with groups of size 1
# Output larger design matrix with spline bases
# main.groups: groups of main effects (original variables)
# all.groups: include main effects and their copies with interactions
# int.groups: list, the [[g]]'th element contains all group indices that group g appears in

generate_gamsel = function(X, groups, degrees = 3, col.normalize = FALSE) {
  group.labels = unique(groups)
  G = length(group.labels)
  gmax = max(group.labels)
  if (length(degrees) == 1) {
    degrees = rep(degrees, G)
  }
  if (length(degrees) != G) {
    stop("Degrees must be a single number or one for each group")
  }
  X.out = X
  n = nrow(X)
  all.groups = groups
  spline.groups = c()
  linear.groups = list()
  gnew = gmax

  for (g in 1:G) {
    group = unique(groups)[g]
    degree = degrees[g]
    gnew = gnew + 1
    spline.groups[group] = gnew
    linear.groups[[gnew]] = group
    col = groups == group
    if (sum(col) > 1) {
      stop("Currently only supports groups of size 1")
    }
    Xvect = X[, col]
    spline.basis = bs(Xvect, degree = degrees[g], Boundary.knots = c(0,1))
    #X.out = cbind(X.out, X[, col])
    X.out = cbind(X.out, spline.basis)
    all.groups = c(all.groups, rep(gnew, ncol(spline.basis))) #+1))
  }
  
  colnames(X.out) = NULL
  if (col.normalize) {
    X.out = col_normalize(X.out)
  }
  return(list(X=X.out, all.groups=all.groups, spline.groups=spline.groups, linear.groups=linear.groups))
}

# Signal vector generation for gamsel
beta_gamsel = function(groups, all.groups, spline.groups, num.nonzero, num.linear, upper, lower, rand.sign=TRUE, perturb=TRUE) {

  main.group.labels = unique(groups)
  all.group.labels = unique(all.groups)
  num.spline = num.nonzero - num.linear
  if (num.spline <= 0) {
    stop("At least 1 nonzero interaction is required")
  }

  nz.linear = sample(main.group.labels, num.linear, replace=FALSE)
  nz.linear.inds = sapply(all.groups, function(x) is.element(x, nz.linear))
  remaining.linear = setdiff(main.group.labels, nz.linear)
  
  if (length(remaining.linear) < num.spline) {
    stop("Not enough variables to form desired number of spline groups")
  }
  
  spline.mains = sample(remaining.linear, num.spline, replace=FALSE)
  nz.spline = c()
  nz.spline.inds = rep(FALSE, length(all.groups))

  for (i in 1:num.spline) {
    this.spline = spline.groups[spline.mains[i]]
    nz.spline = c(nz.spline, this.spline)
    this.spline.inds = which(this.spline == all.groups)
    nz.spline.inds[this.spline.inds] = TRUE
  }

  all.active = c(nz.linear, spline.mains, nz.spline)
  nz.groups = sort(c(nz.linear, nz.spline))
  nz.inds = as.logical(nz.linear.inds + nz.spline.inds)
  if (length(nz.groups) != num.nonzero) {
    stop("Incorrent number of nonzero groups")
  }
  
  magnitudes = sample(seq(from=upper, to=lower, length=num.nonzero))
  beta = rep(0, length(all.groups))

  for (i in 1:num.nonzero) {
    ind.i = all.groups == nz.groups[i]
    beta[ind.i] = magnitudes[i]
  }
  
  if (rand.sign) {
    signs = sample(c(-1, 1), sum(nz.inds), replace = TRUE)
    beta[nz.inds] = signs*beta[nz.inds]
  }

  # Normalize coeff across group for fair comparison with non-grouped vars
  # Normalize interaction terms separately
  
  for (g in nz.groups) {
    group = g == all.groups
    # Linear group, normalize by size
    # (For now, linear groups only have size 1 anyway...)
    if (is.element(g, main.group.labels)) {
      gs = sum(group)
      beta[group] = beta[group]/sqrt(gs)
      bg.norm = sqrt(sum(beta[group]^2))
      
      if (perturb) {
        beta[group] = beta[group] + rnorm(gs)/2
        bg.new.norm = sqrt(sum(beta[group]^2))
        beta[group] = beta[group] * bg.norm / bg.new.norm
      }
      
    } else {
      # Assume first column is linear part
##       group = nz.inds & group # Not necessary?
##       linear.part = nz.linear.inds & group # This will be all 0
##       spline.part = nz.spline.inds & group
##       s.main = sum(main.part)
##       s.int = sum(int.part)
##       beta[main.part] = beta[main.part]/sqrt(s.main)
      linear.part = which(group)[1]
      s.linear = 1
      spline.part = which(group)[-1]
      s.spline = length(spline.part)
      beta[spline.part] = sqrt(2)*beta[spline.part]/sqrt(length(spline.part))
      bg.linear.norm = sqrt(sum(beta[linear.part]^2))
      bg.spline.norm = sqrt(sum(beta[spline.part]^2))
      
      if (perturb) {
        beta[linear.part] = beta[linear.part] + rnorm(s.linear)/2
        bg.main.new.norm = sqrt(sum(beta[linear.part]^2))
        beta[linear.part] = beta[linear.part] * bg.linear.norm / bg.main.new.norm
        
        beta[spline.part] = beta[spline.part] + rnorm(s.spline)/2
        bg.spline.new.norm = sqrt(sum(beta[spline.part]^2))
        beta[spline.part] = beta[spline.part] * bg.spline.norm / bg.spline.new.norm
      }
      
    }
  }

  return(list(beta=beta, true.splines=nz.spline, true.linear = nz.linear, true.active = nz.groups, all.active = all.active))
}
