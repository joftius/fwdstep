
#########     ##       ##       ##
    #       #   #      #  #   #   #
    #       #   #      #  #   #   #
    #        # #       ##      # #

## (1) Fix estimation of beta (keep original X)
## (2) Compare to something natural
## (4) Categorical data?

### Simulation ideas
## Power/performance as function of SNR
## Group sizes differ
# pvals, stepwise with signal in larger vs smaller groups, etc

library(MASS)      # for function 'ginv'
library(dichromat) # for friendly color schemes
library(SGL)       # for comparison

## Functions

####################
# To generate data #
####################

# Fixed group sizes, gaussian design
simulate_fixed = function(n, g, k, orthonormal=TRUE, beta=0) {
  p = g*k
  X = matrix(rnorm(n*p), nrow=n)
  groups = list()
  for (group in 1:g) {
    groups[[group]] = k*(group-1) + 1:k
  }
  if (length(beta) == 1) {
    beta = matrix(0*1:p, ncol=1)
  }
  if (orthonormal == TRUE) {
    for (group in groups) {
      X[,group] = svd(X[,group])$u
    }
  }
  Y = X %*% beta + rnorm(n)
  weights = rep(1,g) * sqrt(k)
  return(list(X=X, Y=Y, groups=groups, weights=weights))
}

# Fixed group sizes, categorical design
simulate_fixed_cat = function(n, g, k, orthonormal=TRUE, beta=0) {
  p = g*k
  X = matrix(nrow=n, ncol=p)
  groups = list()
  for (group in 1:g) {
    groups[[group]] = k*(group-1) + 1:k
    cat.levels = 1
    # Resample until no levels are empty
    while (length(unique(cat.levels)) < k) {
      cat.levels = sample(1:k, n, replace=TRUE)
    }
    cat.binary = unname(model.matrix(~ factor(cat.levels) - 1)[1:n, 1:k])
    if (orthonormal == TRUE) {
      X[1:n, k*(group-1) + 1:k] = cat.binary %*% diag(1/sqrt(colSums(cat.binary)))
    } else {
      X[1:n, k*(group-1) + 1:k] = cat.binary
    }
  }
  if (length(beta) == 1) {
    beta = matrix(0*1:p, ncol=1)
  }
  Y = X %*% beta + rnorm(n)
  weights = rep(1,g) * sqrt(k)
  return(list(X=X, Y=Y, groups=groups, weights=weights))
}


##########################################
# Functions for computing test statistic #
##########################################

# Finds index of group containing var
group_of = function(groups, var) {
  for (i in 1:length(groups)) {
    if (is.element(var[1], groups[[i]])) {
      return(i)
    }
  }
}


# Compute the maximum of the residual process wrt one group
max_resid_proc = function(Xh, wh, Xg, wg, Pg, y) {
  ug = t(Xg) %*% y
  ug = ug / sqrt(sum(ug^2))
  a = t(Xh) %*% (y - Pg %*% y) / wh
  Xug = Xg %*% ug
  b = (wg/wh) * (t(Xh) %*% Xug) / sum(Xug^2)
  ab = sum(a*b)
  norma2 = sum(a^2)
  normb2 = sum(b^2)
  # Rationalized denominator:
  return ((ab + sqrt(norma2*(1-normb2) + ab^2))/(1-normb2))
}


# Compute the test statistic
test_statistic = function(X, Y, groups, weights, active.set=0) {
  inactive.groups = setdiff(1:length(groups), active.set)
  grad = t(X) %*% Y
  Ts = rep(-1,length(groups))
  for (i in inactive.groups) {
    Ts[i] = sqrt(sum(grad[groups[[i]]]^2)) / weights[i]
  }
  imax = which.max(Ts)
  gmax = groups[[imax]]
  X_gmax = X[,gmax]
  #P_gmax = X_gmax %*% ginv(X_gmax)
  P_gmax = X_gmax %*% t(X_gmax)         ### Use orthogonality for now ###
  weight_gmax = weights[imax]
  Ms = c()
  for (i in inactive.groups) {
    if (i != imax) {
      Ms = c(Ms,max_resid_proc(X[,groups[[i]]], weights[i], X_gmax, weight_gmax, P_gmax, Y))
    }
  }
  M = max(Ms)
  u_gmax = grad[gmax]
  u_gmax = u_gmax / sqrt(sum(u_gmax^2))
  f_max = sum(u_gmax * grad[gmax]) / weight_gmax
  var_f_max = sum((X_gmax * u_gmax)^2) / weight_gmax^2
  T = f_max * (f_max - M) / var_f_max
  rank = sum(diag(P_gmax))
  return(list(T=T, lambda1=f_max, lambda2=M, gmax=gmax, imax=imax, weight=weight_gmax, rank=rank))
}


# Avoid numerical overflow
chisq_ratio = function(U, k) {
  numer.p = pchisq(U[1], k, lower.tail=FALSE, log.p=TRUE)
  denom.p = pchisq(U[2], k, lower.tail=FALSE, log.p=TRUE)
  return(exp(numer.p - denom.p))
}


# Compute p-value from test statistic
pval_from_lambdas = function(lambda1, lambda2, weight, rank) {
    U = weight*c(lambda1, lambda2)
    pval = chisq_ratio(U^2, rank)
  if (is.nan(pval)) {
    pval = 1
  }
  return(pval)
}


##################################
# Functions for forward stepwise #
##################################


# Add the next group
add1_lar = function(X, Y, groups, weights, active.set=0, residualize=TRUE, alpha=0.1, origX=X, origY=Y, betahat=0) {
  V = test_statistic(X, Y, groups, weights, active.set)
  V$X = X
  V$Y = Y
  V$betahat = betahat
  V$origX = origX
  V$origY = origY
  V$active.set = active.set
  imax = V$imax
  ginds = groups[[imax]]
  pval = pval_from_lambdas(V$lambda1, V$lambda2, V$weight, V$rank)
  V$pval = pval
  if ((length(betahat) > 1) & (length(active.set) < length(betahat) -1))  {
    V$active.set = setdiff(union(active.set, imax), 0)
    inds = unlist(groups[V$active.set])
    lmnew = lm(origY ~ origX[,inds] - 1)
    betanew = lmnew$coefficients
    V$betahat[inds] = betanew
    V$Y = lmnew$residuals
  }
  if (residualize == TRUE) {
    Xnew = X
    Xg = X[,ginds]
    Hg = diag(1,length(Y)) - Xg %*% t(Xg)
    for (g in setdiff(1:length(groups), imax)) {
      group = groups[[g]]
      Xnew[,group] = svd(Hg %*% X[,group])$u
    }
    V$X = Xnew
  }
  return(V)
}


############################################
# Functions to do a simulation with signal #
############################################


# Calculate statistics like Fdp, etc
selection_results = function(n, g, beta, k, sim.included.all) {
  table.results = matrix(nrow=1, ncol=12)
  sim.included = sim.included.all[[1]]
  q = beta
  fnps = c()
  fdps = c()
  fprs = c()
  fnrs = c()
  for (ll in 1:length(sim.included)) {
    active.set = sim.included[[ll]]
    if (active.set[1] == 0) {
      TP = 0
      FP = 0
      R = 0
    } else {
      TP = sum(active.set <= q)
      FP = sum(active.set > q)
      R = length(active.set)
    }
    
    fdps = c(fdps, ifelse(R == 0, 0, FP/R))  # False discovery proportion
    fprs = c(fprs, FP/(g-q))                 # False positive rate
    fnps = c(fnps, (q-TP)/(g-R))             # False non-discovery proportion
    fnrs = c(fnrs, (q-TP)/q)                 # False negative rate
  }
  table.results[1,] = c(n, g, k, q, mean(fdps), sd(fdps), mean(fprs), sd(fprs),
                        mean(fnps), sd(fnps), mean(fnrs), sd(fnrs))
  colnames(table.results) = c("n", "g", "k", "q", "Fdp", "seFdp", "Fpr", "seFpr", "Fnp", "seFnp", "Fnr", "seFnr")
  return(table.results)
}


# The big simulation function...
sim_forward = function(n=200, g=50, k=5, beta=5, strength=c(1.0, 1.1), maxsteps=20, signaltype="stair",
                       datatype="Gaussian", nsim=100, alpha=0.1, seed=1, echo=TRUE, compare=TRUE) {
  set.seed(seed)
  start.time = Sys.time()
  
  if ((compare == TRUE) & (nsim > 20)) {
    cat(sprintf("Warning: comparing with SGL package will take a long time\n"))
  }
  necho = nsim/10  
  maxdepth = min(c(g-2, floor(n/k), maxsteps))
  p = g*k
  q = k*beta
  
  if (signaltype == "stair") {
    signals = rev(sort(rep(seq(from=max(strength), to=min(strength), length.out = beta), k)))
    betatrue = c(sqrt(2*log(p))*signals, rep(0, p - q))
  } else {
    betatrue = c(rep(sqrt(2*log(p))*max(strength), q), rep(0, p - q))
  }
  
  sim.truebetas = betatrue
  sim.pvals = matrix(1, ncol=maxdepth, nrow=nsim)
  sim.inclmat = sim.pvals*0
  sim.gaps = sim.pvals*0
  sim.betamat = matrix(ncol=p, nrow=nsim)
  sim.included = list()
  sim.trainerrs = c()
  sim.testerrs = c()
  sim.lambdas = matrix(ncol=maxdepth, nrow=nsim)
  sim.mses = c()
  if (compare == TRUE) {
    sgl.betamat = matrix(ncol=p, nrow=nsim)
    sgl.included = list()
    sgl.trainerrs = c()
    sgl.testerrs = c()
    sgl.mses = c()
  }
  t1 = Sys.time()
  for (iter in 1:nsim) {
    if ((iter %% necho == 0) & (echo)) {
      t2 = Sys.time()
      time.left = ((nsim - iter) / necho)*as.numeric(t2-t1)
      cat(sprintf("Sample %4d/%d and time left: %3.1f\n", iter, nsim, time.left))
      t1 = t2
    }
    betahat = rep(0,p)
    if (datatype == "Gaussian") {
      data = simulate_fixed(n,g,k, beta=betatrue)
      testdata = simulate_fixed(n,g,k, beta=betatrue)
    } else if (datatype == "Categorical") {
      data = simulate_fixed_cat(n,g,k, beta=betatrue)
      testdata = simulate_fixed_cat(n,g,k, beta=betatrue)
    }
    Xtest = testdata$X
    Ytest = testdata$Y
    X = data$X
    Y = data$Y
    origY = Y
    origX = X
    groups = data$groups
    weights = data$weights
    active.set = 0
    active.set.stop = 0
    depth = 0
    continue = TRUE
    stopped = FALSE
    while (continue) {
      depth = depth + 1
      V = add1_lar(X, Y, groups, weights, active.set, alpha=alpha, origX=origX, origY=origY, betahat=betahat)
      active.set = V$active.set
      
      ########
      ### old stopping rule
      #######
      if (V$pval >= alpha) {
        stopped = TRUE
      }
      if (stopped == FALSE) {
        active.set.stop = active.set
      }
      betahat = V$betahat
      X = V$X
      Y = V$Y
      sim.pvals[iter,depth] = V$pval
      sim.lambdas[iter,depth] = V$lambda1
      sim.gaps[iter,depth] = V$lambda1 - V$lambda2
      if (depth == maxdepth) {
        continue = FALSE
      }
    }
    
    temp.cand = which(sim.pvals[iter,] < alpha)
    if (length(temp.cand) == 0) {
      sim.included[[iter]] = 0
      active.set = 0
    } else {
      
      print(active.set)
      
      nincl = 1:max(temp.cand)
      active.set = active.set[nincl]
      sim.included[[iter]] = active.set
      sim.inclmat[iter,1:length(active.set)] = active.set
      
      
      print(active.set)
      print("----------------------------------------")
    }
    
    betahat = rep(0,p)
    inds = unlist(groups[active.set])
    betahat[inds] = lm(origY ~ origX[,inds] - 1)$coefficients
    sim.betamat[iter,] = betahat
    sim.mses = c(sim.mses, mean((betatrue-betahat)^2))
    sim.trainerrs = c(sim.trainerrs, mean((origX %*% betahat -origY)^2))
    sim.testerrs = c(sim.testerrs, mean((Ytest - Xtest %*% betahat)^2))
    
    # Compare to cross-validation chosen group lasso from SGL package
    if (compare == TRUE) {
      SGLindex = sort(rep(1:g, k))
      SGLdata = list(x=origX, y=origY)
      SGLfit = cvSGL(data=SGLdata, index=SGLindex, maxit=800, nfold=5, alpha=0)
      lmin = which.min(SGLfit$lldiff)
      lammin = SGLfit$lambdas[lmin]
      SGLbeta = SGLfit$fit
      SGLbeta = SGLbeta$beta[,lmin]
      
      sgl.betamat[iter,] = SGLbeta
      SGLactive = c()
      gcount = 0
      for (group in groups) {
        gcount = gcount + 1
        if(any(SGLbeta[group] != 0)) {
          SGLactive = c(SGLactive, gcount)
        }
      }
      sgl.included[[iter]] = SGLactive
      sgl.trainerrs = c(sgl.trainerrs, mean((origX %*% SGLbeta - origY)^2))
      sgl.testerrs = c(sgl.testerrs, mean((Xtest %*% SGLbeta - Ytest)^2))
      sgl.mses = c(sgl.mses, mean((betatrue-SGLbeta)^2))
    }
    
  }
  sim.pvals.all = list(sim.pvals)
  sim.lambdas.all = list(sim.lambdas)
  sim.betamat.all = list(sim.betamat)
  sim.included.stop = list(sim.included)
  sim.trainerrs.all = list(sim.trainerrs)
  sim.testerrs.all = list(sim.testerrs)
  sim.mses.all = list(sim.mses)
  if (compare == TRUE) {
    sgl.betamat.all = list(sgl.betamat)
    sgl.included.stop = list(sgl.included)
    sgl.trainerrs.all = list(sgl.trainerrs)
    sgl.testerrs.all = list(sgl.testerrs)
    sgl.mses.all = list(sgl.mses)
  }
  # Organize the results
  table.results = matrix(nrow=1, ncol=18)
  table.results[1,1:12] = selection_results(n, g, beta, k, sim.included.stop)
  colnames(table.results) = c("n", "g", "k", "q", "Fdp", "seFdp", "Fpr", "seFpr", "Fnp", "seFnp", "Fnr", "seFnr",
                              "MSE", "seMSE", "TrainErr", "seTrainErr", "TestErr", "seTestErr")
  if (compare == TRUE) {
    sgl.results = matrix(nrow=1, ncol=18)
    sgl.results[1,1:12] = selection_results(n, g, beta, k, sgl.included.stop)
    colnames(sgl.results) = c("n", "g", "k", "q", "Fdp", "seFdp", "Fpr", "seFpr", "Fnp", "seFnp", "Fnr", "seFnr",
                              "MSE", "seMSE", "TrainErr", "seTrainErr", "TestErr", "seTestErr")
  }
  table.results[1,13:18] = c(mean(sim.mses.all[[1]]), sd(sim.mses.all[[1]]),
                             mean(sim.trainerrs.all[[1]]), sd(sim.trainerrs.all[[1]]),
                             mean(sim.testerrs.all[[1]]), sd(sim.testerrs.all[[1]]))
  if (compare == TRUE) {
    sgl.results[1,13:18] = c(mean(sgl.mses.all[[1]]), sd(sgl.mses.all[[1]]),
                             mean(sgl.trainerrs.all[[1]]), sd(sgl.trainerrs.all[[1]]),
                             mean(sgl.testerrs.all[[1]]), sd(sgl.testerrs.all[[1]]))
  }
  
  if (echo) {
    cat(sprintf("Total time: %.2f\n", as.numeric(Sys.time() - start.time)))
  }
  if (compare == TRUE) {
    return(list(step=table.results, sgl=sgl.results, truebetas=sim.truebetas,
                stepbeta=sim.betamat.all, sglbeta=sgl.betamat.all, pvals=sim.pvals.all))
  } else {
    return(list(step=table.results, pvals=sim.pvals, included=sim.inclmat, gaps=sim.gaps))
  }
}

######################################
# Do some simulations and make plots #
######################################

# Marginal plots (jittered)
n=100
g=100
k=5
p=g*k
beta = 10
q=beta*k
strength=c(1.1, 1.05)
nsim = 10

#compare = sim_forward(n, g, k, beta, nsim=nsim, strength=strength, compare=TRUE, maxsteps=40)
compare = sim_forward(n, g, k, beta, nsim=nsim, strength=strength, compare=FALSE, maxsteps=40)




b = beta
incl = compare$included
bb = dim(incl)[2] - b
par(mfrow=c(1,1))
xs = 1:(b+bb)

incl = incl[,xs]
pvm = compare$pvals
pvm = pvm[,xs]
plot(xs, rep(0.2, b+bb), type="l", ylim=c(0,1.05), xlim=c(0.5,b+bb + .5), ylab="p-values", xlab="Order of inclusion",
     main = "g=5/20, n=500, sims=300, stair signal 1 to 1.25")
abline(h=0.2, col="lightgray")
nsims = dim(pvm)[1]
for (j in 1:(b+bb)) {
  pvals = pvm[,j]
  incls = incl[,j]
  cols = ifelse(incls <= b & incls > 0, "black", ifelse(incls <= b, ifelse(j <= b, "black", "blue"), "red"))
  #ifelse(incls <= b & incls > 0, "black", "red")
  cexs = ifelse(incls <= b & j <= b, 0.4, 0.9)
  #pchs = ifelse(incls <= b, ".", 1)
  pchs = sapply(incls, function(g) { if (((g > b) & (j <= b)) | ((g <= b) & (j > b))) return("*") else return(".") })
  pchs = sub("0", ".", as.character(pchs))
  #  if (j > b) {
  #    pchs = "."
  #  }
  points(j + 0.08*rnorm(nsims), pvals, pch=pchs, col=cols, cex=cexs)
}







betatrue = compare$truebetas
vbeta = betatrue
cbeta = betatrue
cbetamat = compare$stepbeta[[1]]
vbetamat = compare$sglbeta[[1]]

cbetamean = colMeans(cbetamat)
cbetasd = apply(cbetamat, 2, sd)

vbetamean = colMeans(vbetamat)
vbetasd = apply(vbetamat, 2, sd)

par(mfrow=c(1,2))
xpts = 1:(4*q)
ymax = max(betatrue) + .5
plot(cbeta[xpts], type="l", ylim=c(-0.1,ymax), xaxt="n", yaxt="n", ylab="Signal", xlab="", main="Forward stepwise")
points(cbetamean[xpts], col="blue")
points(cbetamean[xpts] + cbetasd[xpts], pch=".", col="blue")
points(sapply(cbetamean[xpts] - cbetasd[xpts], function(point) { max(point, 0) }), pch=".", col="blue")

plot(vbeta[xpts], type="l", ylim=c(-0.1,ymax), yaxt="n", xaxt="n", ylab="", xlab="", main="cvSGL")
points(vbetamean[xpts], col="blue")
points(vbetamean[xpts] + vbetasd[xpts], pch=".", col="blue")
points(sapply(vbetamean[xpts] - vbetasd[xpts], function(point) { max(point, 0) }), pch=".", col="blue")
legend("topright", legend=c("True signal", "Mean of est.", "Est. +- MC sd"), col=c("black", "blue", "blue"), pch=c(19, 1 , 46))






par(mfrow=c(1,2))
xpts = 1:(4*q)
ymax = max(cbetamat)
ymin = min(cbetamat)
plot(cbeta[xpts], type="l", lwd=2, ylim=c(ymin,ymax), xaxt="n", yaxt="n", ylab="Signal", xlab="", main="Forward stepwise")
for (row in 1:15) {
  points(cbetamat[row,xpts], pch=".", col="blue")
}

plot(vbeta[xpts], type="l", lwd=2, ylim=c(ymin,ymax), pch=".", yaxt="n", xaxt="n", ylab="", xlab="", main="cvSGL")
for (row in 1:15) {
  points(vbetamat[row,xpts], pch=".", col="blue")
}
legend("topright", legend=c("True signal", "Est. signal"), col=c("black", "blue"), pch=c(19, 46))

