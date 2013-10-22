


pairwise_corrs = function(mat, standardize = FALSE) {
  if (standardize) {
    mat = mat %*% diag(1/sqrt(colSums(mat^2)))
  }
  corrs = c()
  for (i in 1:(ncol(mat)-1)) {
    for (j in (i+1):ncol(mat)) {
      corrs = c(corrs, sum(mat[ , i] * mat[ , j]))
    }
  }
  return(corrs)
}

