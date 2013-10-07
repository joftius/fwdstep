
source('stop_rules.R')

# Input incl.list, p.list, true.active.set
# Use stopping rules stop_first_large, stop_forward, stop_hybrid
# (those input p.list and alpha, output index of last inclusion)

# from fwd_step_test.R results


# input: active.set, true beta, X, Y, X.test, Y.test, etc
# output: prediction error, MSE of beta hat, etc

selection_stats = function(active.set, X, Y, groups, beta, X.test, Y.test, stop.rule, alpha = .1) {
  rule_name = paste("stop_", stop.rule, sep = "")
  if (exists(rule_name, mode = "function")) {
    stop_fun = get(rule_name)
    stop.ind = stop_fun(p.list, alpha)
  } else {
    stop(paste("Misspecified stopping rule:", stop.rule))
  }

  if (stop.ind == 0) {
    stoppped.set = c()
    mse.train = sum(Y.beta^2)
    mse.test = sum(Y.test^2)
    mse.beta = sum(beta^2)
  } else {
    stopped.set = active.set[1:stop.ind]
    groups.active = sapply(groups, function(x) is.element(x, stopped.set))
    fitted.model = lm(Y.beta ~ X[ , groups.active] - 1)
    fitted.beta = fitted.model$coefficients
    mse.train = mean((fitted.model$fitted.values - Y.beta)^2)
    mse.test = mean((predict(fitted.model, newdata = data.frame(X.test)) - Y.test)^2))
    mse.beta = sum((beta[groups.active] - fitted.beta)^2) + sum((beta[setdiff(1:p, which(groups.active))])^2)
  }

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


