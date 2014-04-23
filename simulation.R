source('generate_data.R')
source('fwd_step.R')
source('selection.R')

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
    save=TRUE,
    verbose=FALSE,
    ...) {

    # Load additional files if necessary
    if (type == "glint") {
        # TODO: move glint functions from generate_data.R
        source('glint/generate_glint.R')
    } else if (type == "gamsel") {
        source('gamsel/generate_gamsel.R')
    }

    p = length(groups)
    group.labels = unique(groups)
    g = length(group.labels)
    G = g
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

    P.mat = matrix(nrow=nsim, ncol=max.steps)
    AS.mat = P.mat
    P.mat.b = P.mat
    Chi.mat.b = P.mat
    AS.mat.b = P.mat
    recover.mat = P.mat
    special.recover.mat = P.mat
    ez.recover.mat = P.mat
    ps.recover.mat = P.mat
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
            print("Using fixed design matrix")
            data = fixed.data
            Z = data$X
            # Generate beta here
            
            
        } else if (design == 'gaussian') {
            Z = gaussian_design(n, groups, col.normalize = FALSE, corr = corr)
        } else {
            design_name = paste0(design, "_design")
            if (exists(design_name, mode = "function")) {
                design_fun = get(design_name)
                Z = design_fun(n, groups)
            } else {
                stop(paste("Misspecified design matrix:", design))
            }
        }

        # Construct glint/gamsel design if necessary
        if (type == "default") {
            X = Z
            all.groups = groups
            b.data = beta_staircase(groups, k, upper, lower, cat.groups=cat.groups)
            group.sizes = rle(groups)$lengths
            weights = sqrt(group.sizes)
            
        } else {
            # Nomenclature: "special" means splines or interactions
            # Generate special design matrix
            generate_name = paste0("generate_", type)
            if (exists(generate_name, mode = "function")) {
                generate_fun = get(generate_name)
                data = generate_fun(X=Z, groups=groups, cat.groups=cat.groups, ...)
                X = data$X
                all.groups = data$all.groups
                special.groups = data$special.groups
                default.groups = data$default.groups
                G = length(unique(all.groups))
                # Frobenius normalized, don't need weights
                weights = rep(1, G)
                    
            } else {
                stop(paste("Misspecified simulation type:", type))
            }

            # Attach some functions
            #special_group_of = get(paste0(type, "_special_group_of"))
            default_group_of = get(paste0(type, "_default_group_of"))
            
            # Generate special coefficient vector
            beta_name = paste0("beta_", type)
            beta_fun = get(beta_name)
            b.data = beta_fun(groups, all.groups, special.groups, k=k, num.default=k0, upper=upper, lower=lower, cat.groups=cat.groups)
            # Special active set for glint/gamsel
            true.special = b.data$true.special
                
        } 

        Xnormed = frob_normalize(X, all.groups)
        beta = b.data$beta
        all.active = b.data$all.active
        true.active = b.data$true.active
    
        nz.betas = sapply(unique(all.groups), function(x) sqrt(sum(beta[all.groups==x]^2)))
        nz.betas = nz.betas[nz.betas != 0]
        max.beta = max(nz.betas)
        min.beta = min(nz.betas)
        
        if (ldimS > 1) {
            noise = unwhitener %*% rnorm(n)
        } else if (ldimS == 0) {
            noise = rnorm(n)*sqrt(sigma2)
        } else {
            stop(paste("Misspecified Sigma:", ldimS))
        }
        Y.noiseless = X %*% beta
        mu.max = max(abs(Y.noiseless))
        Y.beta = Y.noiseless + noise
        
#####################################################
        # De-mean before passing to solvers #
        Y.beta = Y.beta - mean(Y.beta)      #
        # Do we really want to do this?     #
#####################################################
        
        # Null results
        results = forward_group(Xnormed, noise, groups=all.groups, weights=weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)
        P.mat[i, ] = results$p.vals
        AS.mat[i, ] = results$active.set

        # Non-null results
        results.b = forward_group(Xnormed, Y.beta, groups=all.groups, weights, Sigma, max.steps = max.steps, cat.groups = cat.groups)

        if (verbose) {
            print(c(upper, lower, mu.max, length(intersect(results.b$active.set[1:k], true.active))/k))
        }
        
        Chi.mat.b[i, ] = results.b$chi.pvals
        P.mat.b[i, ] = results.b$p.vals
        rb.as = results.b$active.set
        recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.active))
        ez.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, all.active))
        AS.mat.b[i, ] = rb.as
        if (type != "default") {
            special.recover.mat[i,] = sapply(rb.as, function(x) is.element(x, true.special))
        }
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
    
    ez.TrueStep = colMeans(ez.recover.mat)
    ez.num.recovered.groups = rowSums(ez.recover.mat[, 1:k])
    TrueStep = colMeans(recover.mat)
    num.recovered.groups = rowSums(recover.mat[, 1:k])
    if (type != "default") {
        num.recovered.special = rowSums(special.recover.mat[, 1:k])
        special.fwd.power = mean(num.recovered.special) / length(true.special)
    }
    fwd.power = mean(num.recovered.groups) / k
    
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
        
        write.csv(out.mats, file = paste0("data/", filename, ".csv"), row.names = FALSE)
        outlist$filename = filename
    }
    return(outlist)
}

