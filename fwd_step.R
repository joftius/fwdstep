
source('group_lasso.R')

update_active_set = function(active.set, group) {

  if (active.set[1] == 0) {
    new.active.set = group
  } else {
    new.active.set = c(active.set, group)
  }
  return(new.active.set)
}

add_group = function(X, Y, groups, weights, active.set = 0) {

  results = group_lasso_knot(X, Y, groups, weights, active.set)
  p.value = pvalue(results$L, results$Mplus, results$Mminus, sqrt(results$var), results$k, sigma=sigma)
  new.active.set = update_active_set(active.set, results$i)

  return(list(test.output = results, p.value = p.value, active.set = new.active.set))
}




forward_group = function(X, Y, groups, weights = 0, max.steps = 0) {

  if (weights == 0) {
    weights = sqrt(rle(groups)$lengths)
  }

  if (max.steps == 0) {
    max.steps = length(unique(groups)) - 1
  }

  active.set = 0
  p.vals = c()
  Ls = c()

  for (i in 1:max.steps) {
    
    output = add_group(X, Y, groups, weights, active.set)
    active.set = output$active.set
    p.vals = c(p.vals, output$p.value)
    Ls = c(Ls, output$test.output[1])
    
  }

  return(list(active.set = active.set, p.vals = p.vals, Ls = Ls))
}
