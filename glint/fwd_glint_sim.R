
source('fwd_step.R')
source('../generate_data.R')
#source('selection.R') not for glinternet
#source('pred_and_estim.R') not for glinternet
#source('fwd_step/coherence.R')

# Note: glint version passes k here instead of beta
fwd_glint_simulation = function(n, Sigma, groups, num.nonzero, lower, upper, nsim,
  max.steps, alpha = .1, design = 'gaussian', corr = 0, categorical = FALSE, predictions = FALSE, coherence = FALSE, plot = FALSE) {

  # Initialize
  p = length(groups)
  g = length(unique(groups))
  G = g*(g+1)/2
  k=num.nonzero
  ldimS = length(dim(Sigma))
  if (ldimS > 1) {
    unwhitener = SigmaSqrt(Sigma)
  } else if (ldimS == 0) {
    Sigma = 1
  }
##   if (g < p) {
##     nz.betas = c()
##     for (g in groups) {
##       bg = sum(beta[groups == g]^2)
##       if (bg > 0) {
##         nz.betas = c(nz.betas, sqrt(bg))
##       }
##     }
##     upper = max(nz.betas)
##     lower = min(nz.betas)
##   } else {
##     upper = max(abs(beta))
##     lower = min(abs(beta[beta != 0]))
##   }
  P.mat = matrix(nrow=nsim, ncol=max.steps)
  AS.mat = P.mat
  P.mat.b = P.mat
  Chi.mat.b = P.mat
  AS.mat.b = P.mat
  recover.mat = P.mat
  int.recover.mat = P.mat
  ez.recover.mat = P.mat
  ps.recover.mat = P.mat
  pred.errs = matrix(0, nrow=3, ncol=3)
  mu.list = c()
  start.time = as.numeric(Sys.time())

  ### Begin main loop ###
  main1 = c()  
  int37 = c()
  
  for (i in 1:nsim) {
    # Monitoring completion time
    if (i %% 10 == 0) {
      elapsed.time = as.numeric(Sys.time()) - start.time
      time.per.iter = elapsed.time/(i-1)
      remaining.time = round((nsim-i)*time.per.iter/60, 3)
      cat(paste("Iteration", i, "of", nsim, "--- time remaining about:", remaining.time, "min\n"))
    }

    if (design == 'gaussian') {
      X = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)
      X.test = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)

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
    all.groups = data$all.groups
    X = data$X
    X.test = data.test$X
#    weights = sqrt(rle(all.groups)$lengths)
#    weights = sqrt(sapply(rle(all.groups)$lengths - 1, function(x) max(1, x)))
    weights = rep(1, G)
    main.groups = data$main.groups
    int.groups = data$int.groups
#source("../generate_data.R")
    beta.data = beta_glinternet(all.groups=all.groups, int.groups=int.groups, num.nonzero=k, upper=upper, lower=lower)
    beta = beta.data$beta
#beta[which(beta!=0)]
#c(sqrt(sum(beta[all.groups==61]^2)),sqrt(sum(beta[all.groups==78]^2)),sqrt(sum(beta[all.groups==94]^2)),sqrt(sum(beta[all.groups==109]^2)))    
    true.active.groups = true_active_groups(all.groups, beta)
    all.active = true.active.groups
    true.ints = beta.data$true.ints
    for (x in true.ints) {
      all.active = c(all.active, main_effects_of(x, int.groups))
    }
    all.active = unname(sort(all.active))
    m=k/3
    ptl.set = unique(all.groups)
    r = length(ptl.set)
    
    # Construct response
    if (ldimS <= 1) {
      Y = rnorm(n)*sqrt(Sigma)
      Y.t = rnorm(n)*sqrt(Sigma)
    } else {
      Y = unwhitener %*% rnorm(n)
      Y.t = unwhitener %*% rnorm(n)
    }

    Y.noiseless = X %*% beta
    mu.max = max(abs(Y.noiseless))
    mu.lc = sum(abs(Y.noiseless)>upper)
    Y.noiseless.test = X.test %*% beta
    Y.beta = Y.noiseless + Y
    Y.test = Y.noiseless.test + Y.t
    X = col_normalize(X)
    X.test = col_normalize(X.test)
    X = frob_normalize(X, all.groups)
    X.test = frob_normalize(X.test, all.groups)    
    main1 = c(main1, abs(t(X[,1]) %*% Y.beta))
    z=t(X[,all.groups==61])%*%Y.beta
    int37 = c(int37, sqrt(sum(z^2)))

    # Null results
    results = forward_group(X, Y, groups=all.groups, weights=weights, Sigma, max.steps = max.steps)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    # Non-null results
    results.b = forward_group(X, Y.beta, groups=all.groups, weights, Sigma, max.steps = max.steps)
print(c(upper, lower, mu.max, mu.lc, length(intersect(results.b$active.set[1:num.nonzero], true.active.groups))/num.nonzero))
    
    Chi.mat.b[i, ] = results.b$chi.pvals
    P.mat.b[i, ] = results.b$p.vals
    rb.as = results.b$active.set
    recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.active.groups))
    int.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.ints))
    ez.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, all.active))

    AS.mat.b[i, ] = rb.as
    already.counted = c()
    cg = 1
    for (ag in rb.as) {
      if (ag <= p) {
        already.counted = union(already.counted, ag)
      } else {
        already.counted = union(already.counted, c(main_effects_of(ag, int.groups), ag))
      }
      ps.recover.mat[i, cg] = length(intersect(already.counted, all.active))/max(length(all.active),1)
      cg = cg + 1
    }

  }

  ez.TrueStep = colMeans(ez.recover.mat)
  ez.num.recovered.groups = rowSums(ez.recover.mat[, 1:num.nonzero])
  TrueStep = colMeans(recover.mat)
  num.recovered.groups = rowSums(recover.mat[, 1:num.nonzero])
  num.recovered.ints = rowSums(int.recover.mat[, 1:num.nonzero])
  fwd.power = mean(num.recovered.groups) / num.nonzero
  int.fwd.power = mean(num.recovered.ints) / length(true.ints)
  ez.fwd.power = mean(ez.num.recovered.groups) / num.nonzero

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

    Chivals = colMeans(Chi.mat.b)
    Chivals.bar <- apply(Chi.mat.b, 2, function(col) quantile(col, probs = bar.quantiles))
    Chivals.point <- apply(Chi.mat.b, 2, function(col) quantile(col, probs = point.quantiles))
    
    xax <- 1:max.steps - 0.15
    nxax <- xax + 0.3
    plot.main <- paste0("n: ", n, ", g: ", g, ", Signal: ", lower.coeff, "/", upper.coeff,  ", Power: ", round(fwd.power, 3), ", Int. Power: ", round(int.fwd.power, 3))

    plot(xax, TrueStep, type = "l", main = plot.main, xlab = "Step", ylab = "", ylim = c(-.05,1.05), xlim = c(min(xax) - .2, max(nxax) + .2), lwd=2)
    points(xax, ez.TrueStep, type = "l", main = plot.main, xlab = "Step", ylab = "", ylim = c(-.05,1.05), xlim = c(min(xax) - .2, max(nxax) + .2), lwd=2, lty=3)
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

    points(xax, Chivals, col="green")
    arrows(xax, Chivals.bar[1, ], xax, Chivals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "green")
    points(xax, Chivals.point[1, ], col = "green", pch = 24, cex = .5)
    points(xax, Chivals.point[2, ], col = "green", pch = 25, cex = .5)

  }

  outlist = list(null.p = P.mat, signal.p = P.mat.b, active.set = AS.mat.b, true.step = recover.mat, psr.mat = ps.recover.mat, m1 = num.nonzero, fwd.power = fwd.power, main1=main1, int37=int37)

  if (predictions) {
    outlist[["pred.errs"]] = pred.errs
  }

  if (coherence) {
    outlist[["coherence"]] = mu.list
  }
  
  return(outlist)
}


print(warnings())

