
source('stop_rules.R')

# Input incl.list, p.list, true.active.set
# Use stopping rules stop_first_large, stop_forward, stop_hybrid
# (those input p.list and alpha, output index of last inclusion)

# from fwd_step_test.R results


# input: active.set, true beta, X, Y, X.test, Y.test, etc
# output: prediction error, MSE of beta hat, etc

pred_est_stats = function(p.list, active.set, X, Y, groups, beta, X.test, Y.test, stop.rule, alpha = .1) {
  p = length(groups)
  rule_name = paste("stop_", stop.rule, sep = "")
  if (exists(rule_name, mode = "function")) {
    stop_fun = get(rule_name)
    stop.ind = stop_fun(p.list, alpha)
  } else {
    stop(paste("Misspecified stopping rule:", stop.rule))
  }
  if (stop.ind == 0) {
    stoppped.set = c()
    mse.train = sum(Y^2)
    mse.test = sum(Y.test^2)
    mse.beta = sum(beta^2)
  } else {
    stopped.set = active.set[1:stop.ind]
    groups.active = sapply(groups, function(x) is.element(x, stopped.set))
    fitted.model = lm(Y ~ X[ , groups.active] - 1)
    fitted.beta = fitted.model$coefficients
    mse.train = mean((fitted.model$fitted.values - Y)^2)
    mse.test = mean((predict(fitted.model, newdata = data.frame(X.test[ , groups.active])) - Y.test)^2)
    mse.beta = sum((beta[groups.active] - fitted.beta)^2) + sum(beta[setdiff(1:p, which(groups.active))]^2)
  }
  return(list(mse.train = mse.train, mse.test = mse.test, mse.beta = mse.beta))
}


# apply above function with all stopping rules
sim_pred_est_stats = function(p.list, active.set, X, Y, groups, beta, X.test, Y.test, alpha) {
  stop.rules = c("first", "forward", "hybrid")
  output = matrix(Inf, nrow = 3, ncol = 1)
  for (stop.rule in stop.rules) {
    output = cbind(output, pred_est_stats(p.list, active.set, X, Y, groups, beta, X.test, Y.test, stop.rule, alpha))
  }
  output = output[ , 1 + 1:length(stop.rules)]
  colnames(output) = stop.rules
  return(output)
}
