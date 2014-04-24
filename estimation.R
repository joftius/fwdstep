
source('stop_rules.R')

# Input incl.list, p.list, true.active.set
# Use stopping rules stop_first_large, stop_forward, stop_hybrid
# (those input p.list and alpha, output index of last inclusion)

# from fwd_step_test.R results


# input: active.set, true beta, X, Y, X.test, Y.test, etc
# output: prediction error, MSE of beta hat, etc

est_stat = function(p.list, active.set, X, Y, groups, beta, X.test, Y.test, stop.rule, alpha = .1) {

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
  
  return(c(mse.train, mse.test, mse.beta))
}


# apply above function with all stopping rules
estimation_stats = function(p.list, active.set, X, Y, groups, beta, X.test, Y.test, alpha = .1) {
  stop.rules = c("first", "forward", "last")
  output = matrix(Inf, nrow = 1, ncol = 3)
  for (stop.rule in stop.rules) {
    rule.err = est_stat(p.list, active.set, X, Y, groups, beta, X.test, Y.test, stop.rule, alpha)
    output = rbind(output, rule.err)
    #print(rule.err)
  }
  output = output[-1, ]
  rownames(output) = stop.rules
  colnames(output) = c("RSS", "Test Error", "textMSE hat beta")
  return(output)
}
