# Tests

source("generate_data.R")

main = function() {
  #n = 20
  #p = 10
  groups = c(1,1,2,2,2,3,3,4,4,5,6,6,6,6,7,8)
  #weights = c(2,2.5,2,2,1.4)
  num.nonzero = 4
  upper = 6
  lower = 5
  beta = beta_staircase(groups, num.nonzero, upper, lower)
  plot(beta, main = "Default")
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.within=TRUE)
  plot(beta, main = "Groups shifted")
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE)
  plot(beta, main = "Random sign")
  beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign=TRUE, rand.within=TRUE)
  plot(beta, main = "Shifted and random sign")
  beta = beta_staircase(groups, num.nonzero, upper, lower, permute=TRUE)
  plot(beta, main = "Permuted")
  beta = beta_staircase(groups, num.nonzero, upper, lower, perturb=TRUE)
  plot(beta, main = "Perturbed")

  n = 100
  groups = sort(c(groups, 5, 7, rep(8, 10)))

  print(groups)
  
  for (test.iter in 1:100) {
    X = categorical_design(n, groups)
    if (check_categorical(X, groups) != TRUE) {
      stop("Categorical design error")
    }
  }
  for (test.iter in 1:20) {
    X = gaussian_design(n, groups, ortho.within = TRUE)
    X.c = categorical_design(n, groups, ortho.within = TRUE)
    if (check_ortho(X, groups) != TRUE) {
      stop("Non-orthogonal gaussian design matrix")
    }
    if (check_ortho(X.c, groups) != TRUE) {
      stop("Non-orthogonal categorical design matrix")
    }
  }
}

check_categorical = function(X, groups) {
  for (g in 1:max(groups)) {
    group = groups == g
    row.sums = rowSums(X[ , group])
    if (any(row.sums != 1)) {
      stop(paste("In group", g, "failure in row(s)", which(row.sums != 1)))
    }
  }
  return(TRUE)
}

check_ortho = function(X, groups, tol = 1e-08) {
  for (g in 1:max(groups)) {
    group = groups == g
    XtX = t(X[ , group]) %*% X[ , group]
    if (max(abs(XtX - diag(1, sum(group)))) > tol) {
      stop(paste("Non-orthogonal submatrix:", g))
    }
  }
  return(TRUE)
}

pdf("test_generate_data.pdf")

main()

dev.off()

