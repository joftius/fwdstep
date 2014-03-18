
source('fwd_step.R')
source('generate_data.R')
source('selection.R')
source('pred_and_estim.R')
source('fwd_step/coherence.R')


fwd_group_simulation = function(n, Sigma, groups, beta, nsim, max.steps,
  alpha = .1, design = 'gaussian', corr = 0, categorical = FALSE,
  predictions = FALSE, rand.beta = FALSE, coherence = FALSE, plot = FALSE,
  fixed.X=NULL, cat.groups = NULL) {

  # Initialize
  weights = sqrt(rle(groups)$lengths)
  p = length(groups)
  g = length(unique(groups))
  true.active.groups = true_active_groups(groups, beta)
  num.nonzero = length(true.active.groups)
  ldimS = length(dim(Sigma))
  if (ldimS > 1) {
    unwhitener = SigmaSqrt(Sigma)
  } else if (ldimS == 0) {
    Sigma = 1
  }   

  if (g < p) {
    nz.betas = c()
    for (g in groups) {
      bg = sum(beta[groups == g]^2)
      if (bg > 0) {
        nz.betas = c(nz.betas, sqrt(bg))
      }
    }
    upper = max(nz.betas)
    lower = min(nz.betas)
  } else {
    upper = max(abs(beta))
    lower = min(abs(beta[beta != 0]))
  }
  P.mat = matrix(nrow=nsim, ncol=max.steps)
  AS.mat = P.mat
  P.mat.b = P.mat
  Chi.mat.b = P.mat
  AS.mat.b = P.mat
  recover.mat = P.mat
  ps.recover.mat = P.mat
  pred.errs = matrix(0, nrow=3, ncol=3)
  mu.list = c()
  start.time = as.numeric(Sys.time())

## track signal reconstruction?  
## don't do this for now
#  if (rand.beta == FALSE) {
#    # track fitted values
#    beta.mat = matrix(0, nrow=nsim, ncol=p)
#  }

  ### Begin main loop ###
  for (i in 1:nsim) {
    # Monitoring completion time
    if (i %% 10 == 0) {
      elapsed.time = as.numeric(Sys.time()) - start.time
      time.per.iter = elapsed.time/(i-1)
      remaining.time = round((nsim-i)*time.per.iter/60, 3)
      cat(paste("Iteration", i, "of", nsim, "--- time remaining about:", remaining.time, "min\n"))
    }

    # Randomized signal
    if (rand.beta == TRUE) {
      beta = beta_staircase(groups, num.nonzero, upper, lower,
        rand.sign = TRUE, permute = TRUE, perturb = TRUE)
      true.active.groups = true_active_groups(groups, beta)
    }

    # Construct design matrix
    if (length(dim(fixed.X)) == 2) {
      X = fixed.X
      # Do not use predictions yet-- need to split to training/test
      X.test = X
    } else {
      if (categorical == TRUE) {
        X = categorical_design(n, groups, ortho.within = FALSE)
        X.test = categorical_design(n, groups, ortho.within = FALSE)
      } else {
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
      }
    }

    # Construct response
    if (ldimS <= 1) {
      Y = rnorm(n)*sqrt(Sigma)
      Y.t = rnorm(n)*sqrt(Sigma)
    } else {
      Y = unwhitener %*% rnorm(n)
      Y.t = unwhitener %*% rnorm(n)
    }
    Y.noiseless = X %*% beta
    Y.noiseless.test = X.test %*% beta
    Y.beta = Y.noiseless + Y
    Y.test = Y.noiseless.test + Y.t
#    X = frob_normalize(X, groups)
#    X.test = frob_normalize(X.test, groups)
#    Y = Y - mean(Y)
#    Y.beta = Y.beta - mean(Y.beta)
    X = col_normalize(X)
    X.test = col_normalize(X.test)

    # Null results
    results = forward_group(X, Y, groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    # Non-null results
    results.b = forward_group(X, Y.beta, groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
print(c(upper, lower, length(intersect(results.b$active.set[1:num.nonzero], true.active.groups))/num.nonzero))
    
    Chi.mat.b[i, ] = results.b$chi.pvals
    P.mat.b[i, ] = results.b$p.vals
    rb.as = results.b$active.set
    AS.mat.b[i, ] = results.b$active.set
    recover.mat[i, ] = sapply(results.b$active.set, function(x)
                 is.element(x, true.active.groups))
    for (cg in 1:length(rb.as)) {
      ps.recover.mat[i, cg] = length(intersect(rb.as[1:cg], true.active.groups))/num.nonzero
    }

    # Tracking prediction error and other things
    # stop.rules must agree with those in sim_pred_est_stats
    if (predictions) {
      if (categorical) {
        stop("Tracking prediction error for categorical design is unsupported")
      }
      stop.rules = c("first", "forward", "last")
      pred.errs = pred.errs + sim_pred_est_stats(results.b$p.vals,
        results.b$active.set, X, Y.beta, groups, beta, X.test, Y.test, alpha)
    }

    # Compute coherence of X
    if (coherence) {
      correlations = pairwise_corrs(X)
      mu = max(abs(correlations))
      mu.list = c(mu.list, mu)
    }

## signal reconstruction    
#    if (rand.beta == FALSE) {
#      #beta.mat[i, groups.active] = fitted.beta
#    }
  }

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

    Chivals = colMeans(Chi.mat.b)
    Chivals.bar <- apply(Chi.mat.b, 2, function(col) quantile(col, probs = bar.quantiles))
    Chivals.point <- apply(Chi.mat.b, 2, function(col) quantile(col, probs = point.quantiles))

    
    xax <- 1:max.steps - 0.15
    nxax <- xax + 0.2
    cxax = xax + 0.4
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

    points(cxax, Chivals, col="green")
    arrows(cxax, Chivals.bar[1, ], cxax, Chivals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "green")
    points(cxax, Chivals.point[1, ], col = "green", pch = 24, cex = .5)
    points(cxax, Chivals.point[2, ], col = "green", pch = 25, cex = .5)

  }

  outlist = list(null.p = P.mat, signal.p = P.mat.b, active.set = AS.mat.b, true.step = recover.mat, psr.mat = ps.recover.mat, m1 = num.nonzero, fwd.power = fwd.power)

  if (predictions) {
    outlist[["pred.errs"]] = pred.errs
  }

  if (coherence) {
    outlist[["coherence"]] = mu.list
  }
  
  return(outlist)
}

print(warnings())

