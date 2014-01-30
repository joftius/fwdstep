
source('fwd_step.R')
source('../generate_data.R')
#source('selection.R') not for glinternet
#source('pred_and_estim.R') not for glinternet
#source('fwd_step/coherence.R')


# Note: glint version passes k here instead of beta
fwd_glint_simulation = function(n, sigma, groups, num.nonzero, lower, upper, nsim,
  max.steps, alpha = .1, design = 'gaussian', corr = 0, categorical = FALSE, predictions = FALSE, coherence = FALSE, plot = FALSE) {

  # Initialize
  weights = sqrt(rle(groups)$lengths)
  p = length(groups)
  g = length(unique(groups))
  k=num.nonzero
  P.mat = matrix(nrow=nsim, ncol=max.steps)
  AS.mat = P.mat
  P.mat.b = P.mat
  AS.mat.b = P.mat
  recover.mat = P.mat
  pred.errs = matrix(0, nrow=3, ncol=3)
  mu.list = c()
  start.time = as.numeric(Sys.time())

  ### Begin main loop ###
  for (i in 1:nsim) {

    # Monitoring completion time
    if (i %% 10 == 0) {
      elapsed.time = as.numeric(Sys.time()) - start.time
      time.per.iter = elapsed.time/(i-1)
      remaining.time = round((nsim-i)*time.per.iter/60, 3)
      cat(paste("Iteration", i, "of", nsim, "--- time remaining about:", remaining.time, "min\n"))
    }

    if (design == 'gaussian') {
      X = gaussian_design(n, groups, corr = corr)
      X.test = gaussian_design(n, groups, corr = corr)

    } else {
      design_name = paste0(design, "_design")
      if (exists(design_name, mode = "function")) {
        design_fun = get(design_name)
        X = design_fun(n, groups)
        X.test = design_fun(n, groups)
      } else {
        stop(paste("Misspecified design matrix:", design))
      }
    }
    # Glinternet specific part
    data = generate_glinternet(X, groups)
    data.test = generate_glinternet(X.test, groups)
    X = data$X
    X.test = data.test$X
    all.groups = data$all.groups
    main.groups = data$main.groups
    int.groups = data$int.groups
    beta.data = beta_glinternet(all.groups=all.groups, int.groups=int.groups, num.nonzero=k, upper=upper, lower=lower)
    beta = beta.data$beta
    true.active.groups = true_active_groups(all.groups, beta)
    true.ints = beta.data$true.ints
    m=k/3
    ptl.set = unique(all.groups)
    r = length(ptl.set)
    true.groups = c(ptl.set[1:m], ptl.set[(r-2*m+1):r])
    
    # Construct response
    Y = rnorm(n)*sigma
    Y.noiseless = X %*% beta
    Y.noiseless.test = X.test %*% beta
    Y.beta = Y.noiseless + Y
    Y.test = Y.noiseless.test + rnorm(n)*sigma

    # Null results
    results = forward_group(X, Y, groups=all.groups, weights, sigma, max.steps = max.steps)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    # Non-null results
    results.b = forward_group(X, Y.beta, groups=all.groups, weights, sigma, max.steps = max.steps)
    P.mat.b[i, ] = results.b$p.vals
    AS.mat.b[i, ] = results.b$active.set
    recover.mat[i, ] = sapply(results.b$active.set, function(x)
                 is.element(x, true.active.groups))

  }


# HERE
  
  TrueStep = colMeans(recover.mat)
  num.recovered.groups = rowSums(recover.mat[, 1:num.nonzero])
  fwd.power = mean(num.recovered.groups) / num.nonzero

  # Convert sum to mean
  pred.errs = pred.errs / nsim

  # Plotting
  if (plot == TRUE) {
    bar.quantiles <- c(.25, .75)
    point.quantiles <- c(.05, .95)
    null.Pvals <- colMeans(P.mat)
    null.Pvals.bar <- apply(P.mat, 2, function(col) quantile(col, probs = bar.quantiles))
    null.Pvals.point <- apply(P.mat, 2, function(col) quantile(col, probs = point.quantiles))

    Pvals = colMeans(P.mat.b)
    Pvals.bar <- apply(P.mat.b, 2, function(col) quantile(col, probs = bar.quantiles))
    Pvals.point <- apply(P.mat.b, 2, function(col) quantile(col, probs = point.quantiles))
    
    xax <- 1:max.steps - 0.15
    nxax <- xax + 0.3
    plot.main <- paste0("n: ", n, ", g: ", g, ", Signal: ", lower.coeff, "/", upper.coeff,  ", k-Oracle power: ", round(fwd.power, 3))

    plot(xax, TrueStep, type = "l", main = plot.main, xlab = "Step", ylab = "", ylim = c(-.05,1.05), xlim = c(min(xax) - .2, max(nxax) + .2), lwd=2)
    abline(v = num.nonzero, lty = "dotted")
    
    points(nxax, null.Pvals, col="orangered")
    arrows(nxax, null.Pvals.bar[1, ], nxax, null.Pvals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "orangered")
    points(nxax, null.Pvals.point[1, ], col = "orangered", pch = 24, cex = .5)
    points(nxax, null.Pvals.point[2, ], col = "orangered", pch = 25, cex = .5)
  
    points(xax, Pvals, col="blue")
    arrows(xax, Pvals.bar[1, ], xax, Pvals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "blue")
    points(xax, Pvals.point[1, ], col = "blue", pch = 24, cex = .5)
    points(xax, Pvals.point[2, ], col = "blue", pch = 25, cex = .5)

  }

  outlist = list(null.p = P.mat, signal.p = P.mat.b, active.set = AS.mat.b, true.step = recover.mat, m1 = num.nonzero, fwd.power = fwd.power)

  if (predictions) {
    outlist[["pred.errs"]] = pred.errs
  }

  if (coherence) {
    outlist[["coherence"]] = mu.list
  }
  
  return(outlist)
}

print(warnings())

