

source('gamsel/fwd_step.R')
source('gamsel/generate_gamsel.R')
source('generate_data.R')
#source('selection.R') not for glinternet
#source('pred_and_estim.R') not for glinternet
#source('fwd_step/coherence.R')

fwd_gamsel_simulation = function(n, Sigma, groups, num.nonzero, num.linear, lower, upper, nsim,  max.steps, design = 'gaussian', corr = 0, categorical = FALSE, predictions = FALSE, coherence = FALSE, fixed.data = NULL, cat.groups = NULL, alpha = .1) {

  # Initialize
  p = length(groups)
  group.labels = unique(groups)
  g = length(group.labels)
  ugsizes = sort(unique(rle(groups)$lengths))
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
  spline.recover.mat = P.mat
  ez.recover.mat = P.mat
  ps.recover.mat = P.mat
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

    # Construct response
    if (length(names(fixed.data)) > 0) {
      print("Using fixed design matrix")
      data = fixed.data
      # Do not use predictions yet-- need to split to training/test
#      data.test = fixed.data
    } else {
      if (design == 'gaussian') {
        X = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)
#        X.test = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)
      } else {
        design_name = paste0(design, "_design")
        if (exists(design_name, mode = "function")) {
          design_fun = get(design_name)
          X = design_fun(n, groups)
#          X.test = design_fun(n, groups)
        } else {
          stop(paste("Misspecified design matrix:", design))
        }
      }
      data = generate_gamsel(X, groups)
#      data.test = generate_glinternet(X.test, groups)    
    }
    X = data$X
#    X.test = data.test$X
    all.groups = data$all.groups
    G = length(unique(all.groups))
#    weights = sqrt(rle(all.groups)$lengths)
#    weights = sqrt(sapply(rle(all.groups)$lengths - 1, function(x) max(1, x)))
    weights = rep(1, G)
    spline.groups = data$spline.groups
    linear.groups = data$linear.groups
    beta.data = beta_gamsel(groups, all.groups, spline.groups, num.nonzero=k, num.linear=num.linear, upper, lower)
    beta = beta.data$beta
    all.active = beta.data$all.active
    true.splines = beta.data$true.splines
    nz.betas = sapply(unique(all.groups), function(x) sqrt(sum(beta[all.groups==x]^2)))
    nz.betas = nz.betas[nz.betas != 0]
    max.beta = max(nz.betas)
    min.beta = min(nz.betas)
    true.active.groups = beta.data$true.active
#    true.active.groups = true_active_groups(all.groups, beta)

    if (ldimS <= 1) {
      Y = rnorm(n)*sqrt(Sigma)
#      Y.t = rnorm(n)*sqrt(Sigma)
    } else {
      Y = unwhitener %*% rnorm(n)
#      Y.t = unwhitener %*% rnorm(n)
    }
    Y.noiseless = X %*% beta
    mu.max = max(abs(Y.noiseless))
#    Y.noiseless.test = X.test %*% beta
    Y.beta = Y.noiseless + Y
#    Y.test = Y.noiseless.test + Y.t
    Y.beta = Y.beta - mean(Y.beta)

    #Xnormed = col_normalize(X)
#    X.test = col_normalize(X.test)
    Xnormed = frob_normalize(X, all.groups)
#    X.test = frob_normalize(X.test, all.groups)    

    # Null results
    results = forward_group(Xnormed, Y, groups=all.groups, weights=weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    # Non-null results
    results.b = forward_group(Xnormed, Y.beta, groups=all.groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
print(c(upper, lower, mu.max, length(intersect(results.b$active.set[1:num.nonzero], true.active.groups))/num.nonzero))
    
    Chi.mat.b[i, ] = results.b$chi.pvals
    P.mat.b[i, ] = results.b$p.vals
    rb.as = results.b$active.set
    recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.active.groups))
    spline.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.splines))
    ez.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, all.active))

    AS.mat.b[i, ] = rb.as
    already.counted = c()
    cg = 1
    for (ag in rb.as) {
      if (is.element(ag, group.labels)) {
        already.counted = union(already.counted, ag)
      } else {
        already.counted = union(already.counted, c(ag, linear.groups[[ag]]))
      }
      ps.recover.mat[i, cg] = length(intersect(already.counted, all.active))/max(length(all.active),1)
      cg = cg + 1
    }

  }

  ez.TrueStep = colMeans(ez.recover.mat)
  ez.num.recovered.groups = rowSums(ez.recover.mat[, 1:num.nonzero])
  TrueStep = colMeans(recover.mat)
  num.recovered.groups = rowSums(recover.mat[, 1:num.nonzero])
  num.recovered.splines = rowSums(spline.recover.mat[, 1:num.nonzero])
  fwd.power = mean(num.recovered.groups) / num.nonzero
  spline.fwd.power = mean(num.recovered.splines) / length(true.splines)
#  ez.fwd.power = mean(ez.num.recovered.groups) / num.nonzero

  # Convert sum to mean
  pred.errs = pred.errs / nsim

##   # Plotting
##   if (plot == TRUE) {
##     bar.quantiles <- c(.25, .75)
##     pospline.quantiles <- c(.05, .95)
##     null.Pvals <- colMeans(P.mat)
##     null.Pvals.bar <- apply(P.mat, 2, function(col) quantile(col, probs = bar.quantiles))
##     null.Pvals.point <- apply(P.mat, 2, function(col) quantile(col, probs = pospline.quantiles))

##     Pvals = colMeans(P.mat.b)
##     Pvals.bar <- apply(P.mat.b, 2, function(col) quantile(col, probs = bar.quantiles))
##     Pvals.point <- apply(P.mat.b, 2, function(col) quantile(col, probs = pospline.quantiles))

##     Chivals = colMeans(Chi.mat.b)
##     Chivals.bar <- apply(Chi.mat.b, 2, function(col) quantile(col, probs = bar.quantiles))
##     Chivals.point <- apply(Chi.mat.b, 2, function(col) quantile(col, probs = pospline.quantiles))
    
##     xax <- 1:max.steps - 0.1
##     cxax = xax + 0.2
##     nxax <- xax + 0.4
##     plot.main <- paste0("n: ", n, ", g: ", g, ", Signal: ", lower.coeff, "/", upper.coeff,  ", Power: ", round(fwd.power, 3), ", Spline. Power: ", round(spline.fwd.power, 3))

##     plot(xax, TrueStep, type = "l", main = plot.main, xlab = "Step", ylab = "", ylim = c(-.05,1.05), xlim = c(min(xax) - .2, max(nxax) + .2), lwd=2)
##     points(xax, ez.TrueStep, type = "l", main = plot.main, xlab = "Step", ylab = "", ylim = c(-.05,1.05), xlim = c(min(xax) - .2, max(nxax) + .2), lwd=2, lty=3)
##     abline(v = num.nonzero, lty = "dotted")
    
##     points(nxax, null.Pvals, col="orangered")
##     arrows(nxax, null.Pvals.bar[1, ], nxax, null.Pvals.bar[2, ],
##            code = 3, angle = 90, length = 0, col = "orangered")
##     points(nxax, null.Pvals.point[1, ], col = "orangered", pch = 24, cex = .5)
##     points(nxax, null.Pvals.point[2, ], col = "orangered", pch = 25, cex = .5)
  
##     points(xax, Pvals, col="blue")
##     arrows(xax, Pvals.bar[1, ], xax, Pvals.bar[2, ],
##            code = 3, angle = 90, length = 0, col = "blue")
##     points(xax, Pvals.point[1, ], col = "blue", pch = 24, cex = .5)
##     points(xax, Pvals.point[2, ], col = "blue", pch = 25, cex = .5)

##     points(cxax, Chivals, col="green")
##     arrows(cxax, Chivals.bar[1, ], cxax, Chivals.bar[2, ],
##            code = 3, angle = 90, length = 0, col = "green")
##     points(cxax, Chivals.point[1, ], col = "green", pch = 24, cex = .5)
##     points(cxax, Chivals.point[2, ], col = "green", pch = 25, cex = .5)

##   }

  outlist = list(TrueStep = TrueStep, null.p = P.mat, signal.p = P.mat.b,
    chi.p = Chi.mat.b, active.set = AS.mat.b, true.step = recover.mat,
    psr.mat = ps.recover.mat, m1 = num.nonzero, fwd.power = fwd.power,
    spline.fwd.power = spline.fwd.power, ugsizes = ugsizes,
    max.beta = max.beta, min.beta = min.beta)

  if (predictions) {
    outlist[["pred.errs"]] = pred.errs
  }

  if (coherence) {
    outlist[["coherence"]] = mu.list
  }
  
  return(outlist)
}


print(warnings())

