source("generate_data.R")
source("fwd_step.R")
source("selection.R")

# Main simulation function
# Saves output in data/
run_tpr_simulation = function(
    nsim = 100,
    type = "default",
    design = "gaussian",
    n, groups, kmax,
    k0 = 1, # For glint/gamsel only
    upper, lower,
    Sigma = 1,
    corr = 0,
    fixed.data = NULL,
    cat.groups = NULL,
    alpha = 0.1,
    noisecorr = 0,
    estimation=TRUE,
    save=TRUE,
    plot=TRUE,
    verbose=TRUE,
    ...) {

    #klist = 3*(1:floor(kmax/3))
    klist = c(1,5,10,15,20,25)
    klen = length(klist)

    if (type != "default") {
        # Load additional files if necessary
        if (type == "glint") {
            source("glint/generate_glint.R")
        } else if (type == "gamsel") {
            source("gamsel/generate_gamsel.R")
        }
        # Attach some functions
        #special_group_of = get(paste0(type, "_special_group_of"))
        default_group_of = get(paste0(type, "_default_group_of"))
        beta_name = paste0("beta_", type)
        beta_fun = get(beta_name)
    }
    
    if (estimation) {
#        if (type != "default") {
#            stop("Estimation error in non-default settings not supported")
#        } else {
            source("estimation.R")
#        }
    }

    p = length(groups)
    group.labels = unique(groups)
    g = length(group.labels)
    G = g
    mult = sqrt(2*log(G)/n)
    upper.scaled = upper * mult
    lower.scaled = lower * mult
    ugsizes = sort(unique(rle(groups)$lengths))
    if (g != p) {
        gmaxmin = paste0(c(min(ugsizes), max(ugsizes)), collapse="-")
    }
    
    # Form unwhitening matrix
    if (noisecorr != 0) {
        # overrides any other specified Sigma
        Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
    } 
    ldimS = length(dim(Sigma))
    if (ldimS > 1) {
        unwhitener = SigmaSqrt(Sigma)
    }  else {
        sigma2 = Sigma
        Sigma = diag(rep(sigma2, n))
    }

    # Initialize storage
    tpp.mat = matrix(nrow=nsim, ncol=klen)
    ez.tpp.mat = tpp.mat
    
    bic.tpp = tpp.mat
    bic.steps = tpp.mat
    ric.tpp = tpp.mat
    ric.steps = tpp.mat
    
    estimation.errs = matrix(0, nrow=1, ncol=klen)
    start.time = as.numeric(Sys.time())

    ### Begin main loop ###
    for (i in 1:nsim) {
        # Monitoring completion time
        if ((i %% 10 == 0) && (i < nsim)) {
            elapsed.time = as.numeric(Sys.time()) - start.time
            time.per.iter = elapsed.time/(i-1)
            remaining.time = round((nsim-i)*time.per.iter/60, 3)
            cat(paste("Iteration", i, "of", nsim, "--- time remaining about:", remaining.time, "min\n"))
        }

        # Construct matrix Z of original features
        if (design == "fixed") {
            if (verbose) print("Using fixed design matrix")
            data = fixed.data
            Z = data$X
            # Generate beta here
            
        } else if (design == "gaussian") {
            Z = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)
            if (estimation) {
                Z.test = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)
            }
        } else {
            design_name = paste0(design, "_design")
            if (exists(design_name, mode = "function")) {
                design_fun = get(design_name)
                Z = design_fun(n, groups)
                if (estimation) {
                    Z.test = design_fun(n, groups, col.normalize = FALSE)
                }   
            } else {
                stop(paste("Misspecified design matrix:", design))
            }
        }

        # Construct glint/gamsel design if necessary
        if (type == "default") {

            # TODO: fixed data here?
            X = Z
            if (estimation) {
                X.test = Z.test
            }
            all.groups = groups
            weights = rep(1, G)
            
        } else {
            # Nomenclature: "special" means splines or interactions
            # Generate special design matrix
            if (design == "fixed") {
                data = fixed.data
                
            } else {
                generate_name = paste0("generate_", type)
                if (exists(generate_name, mode = "function")) {
                    generate_fun = get(generate_name)
                    data = generate_fun(X=Z, groups=groups, cat.groups=cat.groups, ...)
                } else {
                    stop(paste("Misspecified simulation type:", type))
                }
            }
            X = data$X
            all.groups = data$all.groups
            special.groups = data$special.groups
            default.groups = data$default.groups
            G = length(unique(all.groups))
            mult = sqrt(2*log(G)/n)
            upper.scaled = upper * mult
            lower.scaled = lower * mult
            if (estimation) {
                data.test = generate_fun(X=Z.test, groups=groups, cat.groups=cat.groups, ...)
                X.test = data.test$X
            }

            # Frobenius normalized, do not need weights
            weights = rep(1, G)
                
        } 

        Xscaled = scale(X, center = TRUE, scale = FALSE)
        Xscaled = frob_normalize(Xscaled, all.groups)
        if (estimation) {
            Xscaled.test = frob_normalize(X.test, all.groups)
        }

        if (ldimS > 1) {
            noise = unwhitener %*% rnorm(n)
            if (estimation) {
                noise.test = unwhitener %*% rnorm(n)
            }
        } else if (ldimS == 0) {
            noise = rnorm(n)*sqrt(sigma2)
            if (estimation) {
                noise.test = rnorm(n)*sqrt(sigma2)
            }
        } else {
            stop(paste("Misspecified Sigma:", ldimS))
        }

        # Loop through sparsity levels
        for (j in 1:klen) {
            
            k = klist[j]
            # Generate beta
            if (type == "default") {
                b.data = beta_staircase(groups, k, upper.scaled, lower.scaled, cat.groups=cat.groups, rand.sign=TRUE, perturb=TRUE)
                group.sizes = rle(groups)$lengths
                #weights = sqrt(group.sizes)
            } else {
                # Generate special coefficient vector
                b.data = beta_fun(groups=groups, all.groups=all.groups, special.groups=special.groups, default.groups = default.groups, k=k, num.default=k0, upper=upper.scaled, lower=lower.scaled, cat.groups=cat.groups)
                # Special active set for glint/gamsel
                true.special = b.data$true.special
            }
            
            beta = b.data$beta
            all.active = b.data$all.active
            true.active = b.data$true.active
    
            nz.betas = sapply(unique(all.groups), function(x) sqrt(sum(beta[all.groups==x]^2)))
            nz.betas = nz.betas[nz.betas != 0]
            max.beta = max(nz.betas)
            min.beta = min(nz.betas)
        
            Y.noiseless = X %*% beta
            mu.max = max(abs(Y.noiseless))
            Y.beta = Y.noiseless + noise
#####################################################
        # De-mean before passing to solvers #
            Y.beta = Y.beta - mean(Y.beta)
        # Do we really want to do this?     #
#####################################################
        
            if (estimation) {
                Y.noiseless.test = X.test %*% beta
                Y.beta.test = Y.noiseless.test + noise.test
                Y.beta.test = Y.beta.test - mean(Y.beta.test)
            }
        
            # Non-null results
            results = forward_group(Xscaled, Y.beta, groups=all.groups, weights, Sigma, max.steps = k, cat.groups = cat.groups, pval=FALSE)
            as = results$active.set

            # Using step() with BIC/RIC
            if (type == "default") {
              if (design == "gaussian") {
                dfmult = log(n)
                steps = min(floor(n/2), max(groups)-1)
                
                bic.as = step.as(X=Xscaled, Y=Y.beta, steps = steps, k=dfmult)
                b.tpp = length(intersect(bic.as, true.active))/k
                bic.tpp[i,j] = b.tpp
                bic.steps[i,j] = length(bic.as)

                dfmult = 2*log(p)
                ric.as = step.as(X=Xscaled, Y=Y.beta, steps = steps, k=dfmult)
                r.tpp = length(intersect(ric.as, true.active))/k
                ric.tpp[i,j] = r.tpp
                ric.steps[i,j] = length(ric.as)

              }
            }
            
            
            if (estimation) {
                # need function here
            }

            # Track results
            tpp = length(intersect(as[1:k], true.active))/k
            tpp.mat[i,j] = tpp
            if (type != "default") {
                ez.tpp = length(intersect(as[1:k], all.active))/k
                ez.tpp.mat[i, j] = ez.tpp
            }

        } # End sparsity loop

        if (verbose) {
            print(round(c(tpp.mat[i,]), 3))
        }
        
    } # End main simulation loop

    # Summarize some results
    TPR = colMeans(tpp.mat)
    bic.TPR = colMeans(bic.tpp)
    bic.size = as.character(apply(bic.steps, 2, median))
    ric.TPR = colMeans(ric.tpp)
    ric.size = as.character(apply(ric.steps, 2, median))

    
    if (type != "default") {
        ez.TPR = colMeans(ez.tpp.mat)
    }
    estimation.errs = estimation.errs / nsim
    
    outlist = list(
        n = n, p = p, g = g, klist = klist,
        TPR = TPR,
        tpp = tpp.mat,
        ugsizes = ugsizes,
        max.beta = max.beta,
        min.beta = min.beta)
    
    if (type != "default") {
        outlist$ez.tpp = ez.tpp.mat
        outlist$ez.TPR = ez.TPR
    }
    if (estimation) {
        outlist$estimation = estimation.errs
    }

    if (save) {
        filename = paste0("tpr_", type, "_", design, "_nsim", nsim,
            "_n", n, "_p", p, "_g", g)
        if (g != p) {
            filename = paste0(filename, "_sizes", gmaxmin)
        }
        filename = paste0(filename, "_k", k)
        if (type != "default") {
            filename = paste0(filename, "_k0", k0)
        }
        filename = paste0(filename, "_lower", lower, "_upper", upper)
        if (corr != 0) {
            filename = paste0(filename, "_corr", corr)
        }
        if (noisecorr != 0) {
            filename = paste0(filename, "_noisecorr", noisecorr)
        }
        
        out.mats = data.frame(iter=1:nsim, tpp=tpp.mat)
        
        filename = gsub(".", "pt", filename, fixed=TRUE)
        write.csv(out.mats, file = paste0("data/", filename, ".csv"), row.names = FALSE)
        outlist$filename = filename
    }

    if (plot) {
      plotfile = paste0("figs/", filename, ".pdf")

      plot.main <- paste0("n:", n, ", p:", p, ", g:", g)
      if (g != p) {
        plot.main = paste0(plot.main, "(", gmaxmin, ")")
      }
      plot.main = paste0(plot.main, 
        ", beta:", upper, "/", lower,
        "(", round(max.beta,1), "/", round(min.beta, 1),
        ")")
      if ((corr != 0) | (noisecorr != 0)) {
        plot.main = paste0(plot.main,
          ", corr:", corr, "/", noisecorr)
      }
  
      pdf(plotfile)
      plot(klist, TPR, type = "l", main = plot.main, xlab = "Sparsity", ylab = "TPP", ylim = c(-0.1, 1.1), lwd = 2, xaxt="n")
      axis(1, at=klist)

      abline(h = c(.9, .7, .5, .3, .1), lty = 2, col = "gray")
      abline(h = c(1, .8, .6, .4, .2, 0), lty = 3, col = "gray")
  
      for (j in 1:klen) {
        xjitter = 0.05*rnorm(nsim)
        yjitter = 0.01*rnorm(nsim)
        points(rep(klist[j], nsim) + xjitter, tpp.mat[, j] + yjitter, pch = ".", col = rgb(red=0, green=0, blue=1), cex = 2)
      }

      points(klist, TPR, type = "l", lwd = 2)
      
      if (type != 'default') {
        points(klist, ez.TPR, type="l", lwd = 2, lty="dashed")
      } else {
        if (design == "gaussian") {
          points(klist, bic.TPR, type="l", lwd = 2, lty="dashed", col="green")
          axis(side=3, at=klist, labels=bic.size, pos=1.02, tck=0,
               lwd=0, col.axis = "green")

          points(klist, ric.TPR, type="l", lwd = 2, lty="dotdash", col="purple")
          axis(side=1, at=klist, labels=ric.size, pos=-0.02, tck=0,
                 lwd=0, col.axis = "purple")
        }
      }

      dev.off()
    }

    print(ric.size)
    print(bic.size)

    return(outlist)
}

