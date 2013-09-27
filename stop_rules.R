
# Input a list of p-values and alpha
# Output largest index for model inclusion (0 if empty model)

stop_first_large = function(p.list, alpha = .1) {

  p.big = which(p.list > alpha)

  if (length(p.big) == 0) {
    return(length(p.list))
  } else {
    k = min(p.big) - 1
    return(k)
  }
}


stop_forward = function(p.list, alpha = .1) {

  k = 0
  sum.log = 0
  log.list = log(1-p.list)
  continue = TRUE
  
  while (continue) {
    
    k = k + 1
    sum.log = sum.log + log.list[k]
    if (k >= length(p.list)) {
      continue = FALSE
    }
    if (sum.log > -k*alpha) {
      continue = FALSE
      k = k - 1
    }
    
    if (k == 0) {
      return(0)
    } else {
      return(k)
    }
  }
}


stop_hybrid = function(p.list, alpha = .1) {
  return(2)
}

