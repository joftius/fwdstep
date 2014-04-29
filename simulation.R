source("generate_data.R")
source("fwd_step.R")
source("selection.R")

# Main simulation function
# Saves output in data/
run_simulation = function(
    nsim = 500,
    type = "default",
    design = "gaussian",
    n, groups, k,
    k0 = 1, # For glint/gamsel only
    upper, lower,
    max.steps = NULL,
    Sigma = 1,
    corr = 0,
    fixed.data = NULL,
    cat.groups = NULL,
    alpha = 0.1,
    noisecorr = 0,
    estimation=FALSE,
    save=TRUE,
    verbose=FALSE,
    ...) {


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
    P.mat = matrix(nrow=nsim, ncol=max.steps)
    AS.mat = P.mat
    P.mat.b = P.mat
    Chi.mat.b = P.mat
    AS.mat.b = P.mat
    recover.mat = P.mat
    special.recover.mat = P.mat
    ez.recover.mat = P.mat
    ps.recover.mat = P.mat
    estimation.errs = matrix(0, nrow=3, ncol=3)
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
            b.data = beta_staircase(groups, k, upper.scaled, lower.scaled, cat.groups=cat.groups, rand.sign=TRUE, perturb=TRUE)
            group.sizes = rle(groups)$lengths
            weights = sqrt(group.sizes)
            
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
            

            
            # Generate special coefficient vector
            b.data = beta_fun(groups=groups, all.groups=all.groups, special.groups=special.groups, default.groups = default.groups, k=k, num.default=k0, upper=upper.scaled, lower=lower.scaled, cat.groups=cat.groups)
            # Special active set for glint/gamsel
            true.special = b.data$true.special
                
        } 

        Xscaled = scale(X, center = TRUE, scale = FALSE)
        Xscaled = frob_normalize(Xscaled, all.groups)
        if (estimation) {
            Xscaled.test = frob_normalize(X.test, all.groups)
        }
        beta = b.data$beta
        all.active = b.data$all.active
        true.active = b.data$true.active
    
        nz.betas = sapply(unique(all.groups), function(x) sqrt(sum(beta[all.groups==x]^2)))
        nz.betas = nz.betas[nz.betas != 0]
        max.beta = max(nz.betas)
        min.beta = min(nz.betas)
        
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
        
        # Null results
        results = forward_group(Xscaled, noise, groups=all.groups, weights=weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
        P.mat[i, ] = results$p.vals
        AS.mat[i, ] = results$active.set

        # Non-null results
        results.b = forward_group(Xscaled, Y.beta, groups=all.groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)

        if (verbose) {
            print(c(upper, lower, mu.max, length(intersect(results.b$active.set[1:k], true.active))/k))
        }

        # Track results
        Chi.mat.b[i, ] = results.b$chi.pvals
        P.mat.b[i, ] = results.b$p.vals
        rb.as = results.b$active.set
        recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.active))
        ez.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, all.active))
        AS.mat.b[i, ] = rb.as

        # Track estimation results
        if (estimation) {
            stop.rules = c("first", "forward", "last")
            estimation.errs = estimation.errs + estimation_stats(
                p.list = results$p.vals,
                active.set = rb.as, X = Xscaled, Y = Y.beta,
                groups = all.groups, beta = beta,
                X.test = Xscaled.test, Y.test = Y.beta.test, alpha = alpha)
        }

        if (type != "default") {
            special.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.special))
        }

        # Track proportion of signal recovered
        already.counted = c()
        cg = 1
        for (ag in rb.as) {
            # Append already.counted
            if (type == "default") {
                already.counted = c(already.counted, ag)
            } else {
                if (is.element(ag, group.labels)) {
                   # Default group was added
                    already.counted = union(already.counted, ag)
                } else {
                    # Special group was added, count default also
                    already.counted = union(already.counted, c(ag, default_group_of(ag, default.groups)))
                }
            }
            # Append proportion of signal recover matrix
            ps.recover.mat[i, cg] = length(intersect(already.counted, all.active))/max(length(all.active),1)
            cg = cg + 1
        }
    }
    # End main simulation loop

    # Summarize some results
    ez.TrueStep = colMeans(ez.recover.mat)
    ez.num.recovered.groups = rowSums(ez.recover.mat[, 1:k])
    TrueStep = colMeans(recover.mat)
    num.recovered.groups = rowSums(recover.mat[, 1:k])
    if (type != "default") {
        num.recovered.special = rowSums(special.recover.mat[, 1:k])
        special.fwd.power = mean(num.recovered.special) / length(true.special)
    }
    fwd.power = mean(num.recovered.groups) / k
    estimation.errs = estimation.errs / nsim
    
    out.mats = data.frame(iter=1:nsim, null = (P.mat), signal = (P.mat.b), chi = (Chi.mat.b), true.step = (recover.mat), prop.signal = (ps.recover.mat))
    outlist = list(
        n = n, p = p, g = g, k = k,
        TrueStep = TrueStep,
        null.p = P.mat,
        signal.p = P.mat.b,
        chi.p = Chi.mat.b,
        active.set = AS.mat.b,
        true.step = recover.mat,
        psr.mat = ps.recover.mat,
        fwd.power = fwd.power,
        ugsizes = ugsizes,
        max.beta = max.beta,
        min.beta = min.beta)
    
    if (type != "default") {
        outlist$special.fwd.power = special.fwd.power
    }
    if (estimation) {
        outlist$estimation = estimation.errs
    }

    if (save) {
        filename = paste0(type, "_", design, "_nsim", nsim,
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

        filename = gsub(".", "pt", filename, fixed=TRUE)
        write.csv(out.mats, file = paste0("data/", filename, ".csv"), row.names = FALSE)
        outlist$filename = filename
    }
    return(outlist)
}

