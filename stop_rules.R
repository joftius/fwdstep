# TODO
# Figure out how to use HybridStop without T_k
# or convert to/calculate T_k from our p-values?

# Input a list of p-values and alpha
# Output largest index for model inclusion (0 if empty model)

stop_first = function(p.list, alpha = .1) {

  p.big = which(p.list > alpha)

  if (length(p.big) == 0) {
    return(length(p.list))
  } else {
    k = min(p.big) - 1
    return(k)
  }
}

# Renyi representation
p_to_q = function(p.list) {
  m = length(p.list)
  y = -log(1-p.list)
  denoms = m + 1 - 1:m
  z = cumsum(y/denoms)
#  Z = rep(0, m)
#  for (i in 1:m) {
#    Z[i] = sum(y[1:i]/(m-1:i+1))
#  }
  return(z)
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
    if (-sum.log/k > alpha) {
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

p_to_exp = function(p.list) {
  
}


stop_hybrid = function(p.list, tau = 5, alpha = .1) {
  q = p_to_q(p.list)
  
  return(2)
}

