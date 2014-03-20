
source('fwd_step.R')
source('generate_data.R')
source('selection.R')
source('pred_and_estim.R')
source('fwd_step/coherence.R')


fwd_group_simulation = function(n, Sigma, groups, beta, nsim, max.steps,
  alpha = .1, design = 'gaussian', corr = 0, categorical = FALSE,
  predictions = FALSE, rand.beta = TRUE, coherence = FALSE,
  fixed.data=NULL, cat.groups = NULL, staircase = TRUE, SNR = 1) {

  # Initialize
  gsizes = rle(groups)$lengths
  ugsizes = sort(unique(gsizes))
  weights = sqrt(gsizes)
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
  beta.list = c()
  start.time = as.numeric(Sys.time())

  R.powdiff = c()
  R.iwin = c()
  R.rwin = c()
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
    if (rand.beta) {
      if (categorical) {
        beta = beta_staircase(groups, num.nonzero, upper, lower,
          permute = TRUE, perturb = TRUE, cat.groups = unique(groups))
      } else {
        beta = beta_staircase(groups, num.nonzero, upper, lower,
          rand.sign = TRUE, permute = TRUE, perturb = TRUE)
      }
      true.active.groups = true_active_groups(groups, beta)
      beta = beta
    }
    # Construct design matrix
    if (length(fixed.data) > 0) {
      print("Using fixed data set")
      X = fixed.data$fixed.X
      X.cat = fixed.data$X.cat
      # Do not use predictions yet-- need to split to training/test
      X.test = X
    } else {
      
      if (categorical == TRUE) {
        catd = categorical_design(n, groups)
        catd.test = categorical_design(n, groups)
        X = catd$X
        X.test = catd.test$X
        X.cat = catd$X.cat
        
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

    # Scaling
    Xnormed = frob_normalize(X, groups)
    #X.test = frob_normalize(X.test, groups)
    
    nz.betas = sapply(unique(groups), function(x) sqrt(sum(beta[groups==x]^2)))
    nz.betas = nz.betas[nz.betas != 0]
    max.beta = max(nz.betas)
    min.beta = min(nz.betas)
    mean.beta = mean(nz.betas)
    beta.list = c(beta.list, mean.beta)

    # True mean formed from original X
    Y.noiseless = X %*% beta
    mu.max = max(abs(Y.noiseless))
    Y.noiseless.test = X.test %*% beta
    Y.beta = Y.noiseless + Y
    Y.test = Y.noiseless.test + Y.t
    # Centering
    Y = Y - mean(Y)
    Y.beta = Y.beta - mean(Y.beta)

    # Null results
    results = forward_group(Xnormed, Y, groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
    P.mat[i, ] = results$p.vals
    AS.mat[i, ] = results$active.set

    # Non-null results
    results.b = forward_group(Xnormed, Y.beta, groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
    rb.as = results.b$active.set
    my.pow = length(intersect(rb.as[1:num.nonzero], true.active.groups))/num.nonzero
    trace = c(upper, lower, mu.max, my.pow)
    
    if (length(cat.groups) > 0) {

      null.fit = lm(Y.beta ~ 0, data=X.cat)
      full.fit = lm(Y.beta ~ . -1, data=X.cat)
      step.fit = step(null.fit, scope=list(lower=null.fit, upper=full.fit), direction = "forward", k = 0, steps = num.nonzero, trace = 0)
      step.groups = gsub("X", "", names(step.fit$coefficients), fixed=TRUE)
      step.groups = sapply(step.groups, function(s) substr(s, 1, nchar(s) - 1))
      step.groups = as.numeric(unique(step.groups))
      Rs.pow = length(intersect(step.groups, true.active.groups))/num.nonzero
      trace = c(trace, Rs.pow)
      R.powdiff = R.powdiff + my.pow - Rs.pow
      R.iwin = c(R.iwin, my.pow > Rs.pow)
      R.rwin = c(R.rwin, my.pow < Rs.pow)
      
    } else {
      
      df = data.frame(Xnormed)
      null.fit = lm(Y.beta ~ 0, data=df)
      full.fit = lm(Y.beta ~ . -1, data=df)
      step.fit = step(null.fit, scope=list(lower=null.fit, upper=full.fit), direction = "forward", k = 0, steps = num.nonzero, trace = 0)
      step.groups = as.numeric(gsub("X", "", names(step.fit$coefficients), fixed=TRUE))
      Rs.pow = length(intersect(step.groups, true.active.groups))/num.nonzero
      trace = c(trace, Rs.pow)
      R.powdiff = R.powdiff + my.pow - Rs.pow
      R.iwin = c(R.iwin, my.pow > Rs.pow)
      R.rwin = c(R.rwin, my.pow < Rs.pow)
      
    }
    
    print(round(trace, 2))
    
    Chi.mat.b[i, ] = results.b$chi.pvals
    P.mat.b[i, ] = results.b$p.vals

    AS.mat.b[i, ] = rb.as
    recover.mat[i, ] = sapply(rb.as, function(x)
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
##       correlations = pairwise_corrs(X)
##       mu = max(abs(correlations))
##       mu.list = c(mu.list, mu)
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
  mean.beta = mean(beta.list)

  outlist = list(TrueStep = TrueStep, null.p = P.mat, signal.p = P.mat.b, chi.p = Chi.mat.b, active.set = AS.mat.b, true.step = recover.mat, psr.mat = ps.recover.mat, m1 = num.nonzero, fwd.power = fwd.power, ugsizes = ugsizes, min.beta = min.beta, max.beta = max.beta, mean.beta)

  if (predictions) {
    outlist[["pred.errs"]] = pred.errs
  }

  if (coherence) {
    outlist[["coherence"]] = mu.list
  }

  print(c(sum(R.iwin), sum(R.rwin), R.powdiff/nsim))
  
  return(outlist)
}

print(warnings())

