
source('stop_rules.R')

# Input incl.list, p.list, true.active.set
# Use stopping rules stop_first_large, stop_forward, stop_last
# (those input p.list and alpha, output index of last inclusion)

# from fwd_step_test.R results


# write a function taking  these inputs and calculating:
# Rob/Alex/Max/Stefan show FDR, V, R, and average power (fraction of non-null rejected)

selection_stats = function(p.list, active.set, true.step, k, stop.rule, alpha = .1) {
  rule_name = paste("stop_", stop.rule, sep = "")
  if (exists(rule_name, mode = "function")) {
    stop_fun = get(rule_name)
    stop.ind = stop_fun(p.list, alpha)
  } else {
    stop(paste("Misspecified stopping rule:", stop.rule))
  }
  if (stop.ind == 0) {
    #stopped.set = c()
    S = 0
  } else {
    #stopped.set = active.set[1:stop.ind]
    S = sum(true.step[1:stop.ind])
  }
  R = stop.ind
  V = R - S
  Fdp = V/max(1, R)
  Power = S/max(1, k)
  return(c(R, Fdp, Power))
}

# names = null.p, signal.p, active.set, true.set
sim_select_stats = function(p.lists, active.sets, true.steps, k, alpha = .1, bic.list = NULL, ric.list = NULL) {
  first.mat = matrix(0, nrow = nrow(p.lists), ncol = 3)
  forward.mat = first.mat
  last.mat = first.mat
  bic.mat = first.mat
  ric.mat = first.mat
  oracle.mat = first.mat
  for (i in 1:nrow(p.lists)) {
    oracle.mat[i, ] = selection_stats(p.lists[i, 1:k], active.sets[i, ], true.steps[i, ], k, stop.rule = "all", alpha)
    first.mat[i, ] = selection_stats(p.lists[i, ], active.sets[i, ], true.steps[i, ], k, stop.rule = "first", alpha)
    forward.mat[i, ] = selection_stats(p.lists[i, ], active.sets[i, ], true.steps[i, ], k, stop.rule = "forward", alpha)
    last.mat[i, ] = selection_stats(p.lists[i, ], active.sets[i, ], true.steps[i, ], k, stop.rule = "last", alpha)
    if (length(bic.list) > 0) {
      bics = bic.list[[i]]
      bic.mat[i, ] = selection_stats(p.list=rep(0, bics[1]), active.sets[i, ], true.step=bics[-1], k, stop.rule = "all", alpha)
    }
    if (length(ric.list) > 0) {
      rics = ric.list[[i]]
      ric.mat[i, ] = selection_stats(p.list=rep(0, rics[1]), active.sets[i, ], true.step=rics[-1], k, stop.rule = "all", alpha)
    }
  }
  oracle = MC_table(oracle.mat)
  first = MC_table(first.mat)
  forward = MC_table(forward.mat)
  last = MC_table(last.mat)
  results = rbind(oracle, first, forward, last)
  if (length(ric.list) > 0) {
    RIC = MC_table(ric.mat)
    results = rbind(results, RIC)
  }
  if (length(bic.list) > 0) {
    BIC = MC_table(bic.mat)
    results = rbind(results, BIC)
  }
  colnames(results) = c("R", "FDP", "TPP")
  rownames(results) = paste(k, rownames(results))
  return(results)
}

MC_table = function(mat) {
  means = round(colMeans(mat), 2)
  means[1] = round(means[1], 1)
  sds = round(apply(mat, 2, sd), 2)
  sds[1] = round(sds[1], 1)
  out = c()
  for (i in 1:length(means)) {
    out = c(out, paste0(means[i], "(", sds[i], ")"))
  }
  return(out)
}
