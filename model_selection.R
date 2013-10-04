
source('stop_rules.R')

# Input incl.list, p.list, true.active.set
# Use stopping rules stop_first_large, stop_forward, stop_hybrid
# (those input p.list and alpha, output index of last inclusion)

# from fwd_step_test.R results


# write a function taking  these inputs and calculating:
# Rob/Alex/Max/Stefan show FDR, V, R, and average power (fraction of non-null rejected)

selection_stats = function(p.list, active.set, true.step, m1, stop.rule, alpha = .1) {
  rule_name = paste("stop_", stop.rule, sep = "")
  if (exists(rule_name, mode = "function")) {
    stop_fun = get(rule_name)
    stop.ind = stop_fun(p.list, alpha)
  } else {
    stop(paste("Misspecified stopping rule:", stop.rule))
  }
#  if (stop.rule == "first") {
#    stop.ind = stop_first(p.list, alpha)
#  } else if (stop.rule == "forward") {
#    stop.ind = stop_forward(p.list, alpha)
#  } else if (stop.rule == "hybrid") {
#    stop.ind = stop_hybrid(p.list, alpha)
  if (stop.ind == 0) {
    stopped.set = c()
    S = 0
  } else {
    stopped.set = active.set[1:stop.ind]
    S = sum(true.step[1:stop.ind])
  }
  R = stop.ind
  V = R - S
  Fdp = V/max(1, R)
  Power = S/max(1, m1)
  return(c(Fdp, R, S, V, Power))
}

# names = null.p, signal.p, active.set, true.set
sim_select_stats = function(p.lists, active.sets, true.steps, m1, alpha = .1) {
  first.mat = matrix(0, nrow = nrow(p.lists), ncol = 5)
  forward.mat = first.mat
  hybrid.mat = first.mat
  for (i in 1:nrow(p.lists)) {
    first.mat[i, ] = selection_stats(p.lists[i, ], active.sets[i, ], true.steps[i, ], m1, stop.rule = "first", alpha)
    forward.mat[i, ] = selection_stats(p.lists[i, ], active.sets[i, ], true.steps[i, ], m1, stop.rule = "forward", alpha)
    hybrid.mat[i, ] = selection_stats(p.lists[i, ], active.sets[i, ], true.steps[i, ], m1, stop.rule = "hybrid", alpha)
  }
  first = colMeans(first.mat)
  forward = colMeans(forward.mat)
  hybrid = colMeans(hybrid.mat)
  results = rbind(first, forward, hybrid)
  colnames(results) = c("Fdp", "R", "S", "V", "Power")
  return(results)
}


