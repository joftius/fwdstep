
source('../generate_data.R')
source('fstep_minimal.R')


##########################
### Modify these lines ###
fname = "p20"           ##
ns = c(200) #, 200, 400)   ##
ps = 40 # floor(ns/2)        ##
nsim = 100               ##
klist = 1:4
#########################
#fname = "p2n"          ##
#ns = c(100) #, 200)    ##
#ps = ceiling(ns*2)     ##
#nsim = 75              ##
#kplus = 10             ##
##########################
corr = 0
##########################

lkl = length(klist)

for (i in 1:length(ns)) {
  mu.mat = matrix(0, nrow=lkl, ncol = nsim)
  pow.mat = mu.mat
  ez.pow.mat = pow.mat
  
  n = ns[i]
  p = ps[i]
  mult = sqrt(2*log(p*(p+1)/2))
  upper.coeff = 2.9
  lower.coeff = 1.5
  upper = upper.coeff*mult
  lower = lower.coeff*mult
  klist = 3*klist
  groups = 1:p

  for (j in 1:lkl) {
    k = klist[j]
    print(k)
    
    for (iter in 1:nsim) {

      data = generate_glinternet(gaussian_design(n, groups, corr = corr), groups)
      X = data$X
      all.groups = data$all.groups
      main.groups = data$main.groups
      int.groups = data$int.groups
      beta.data = beta_glinternet(all.groups=all.groups, int.groups=int.groups, num.nonzero=k, upper=upper, lower=lower)
      beta = beta.data$beta
      true.ints = beta.data$true.ints
      m = k/3
      ptl.set = unique(all.groups)
      r = length(ptl.set)
      true.active.groups = true_active_groups(all.groups, beta)
      Y = X %*% beta + rnorm(n)
      Y = Y - mean(Y)

      # Weights?
      #weights = sqrt(rle(all.groups)$lengths)
      weights = rep(1, max(all.groups)) #c(rep(1, p), rep(sqrt(2), max(all.groups) - p))
      #weights = sqrt(sapply(rle(all.groups)$lengths - 1, function(x) max(1, x)))
      
      # Fit forward stepwise
      active.set = fstep_fit(X=X, Y=Y, groups=all.groups, max.steps=k, weights=weights)

##       c.pow = 0
##       already.counted = c()
##       cg = 1
##       for (ag in active.set) {
##         c.pow = c.pow + true_step_glinternet(ag, p, int.groups, active.set[1:cg], true.active.groups, already.counted)
##         if (ag <= p) {
##           already.counted = union(already.counted, ag)
##         } else {
##           already.counted = union(already.counted, c(main_effects_of(ag, int.groups), ag))
##         }
##         cg = cg + 1
##       }
      
      true.added = intersect(active.set, true.active.groups)
      true.added.ints = intersect(true.added, true.ints)
      false.added = setdiff(active.set, true.active.groups)
      ez.pow.mat[j, iter] = length(true.added.ints)/length(true.ints)
      pow.mat[j, iter] = length(true.added)/k

    }
  }
  
  power.MCavg = rowMeans(pow.mat)
  ez.power.MCavg = rowMeans(ez.pow.mat)

  print(rbind(klist, power.MCavg))

  # Save plot
  plot.main = paste0("n = ", n, ", p = ", p, ", signal strength ", lower.coeff, "/", upper.coeff)
  filename = paste0('../figs/fwdstep/gauss_glint_', fname, "_n", n, "_lower", lower.coeff, "_upper", upper.coeff)
  if (corr != 0) {
    filename = paste0(filename, "_corr", corr)
    plot.main = paste0(plot.main, ", corr = ", corr)
  }
  filename = paste0(filename, ".pdf")
  
  pdf(filename)
  plot(klist, power.MCavg, type = "l", main = plot.main, xlab = "Sparsity", ylab = "Average power", ylim = c(-0.1, 1.1), lwd = 2)
  points(klist, ez.power.MCavg, type = "l", lty = 3, lwd = 2)
  abline(h = c(.9, .7, .5, .3, .1), lty = 2, col = "gray")
  abline(h = c(1, .8, .6, .4, .2, 0), lty = 3, col = "gray")
  
  for (j in 1:lkl) {
    xjitter = 0.05*rnorm(nsim)
    yjitter = 0.01*rnorm(nsim)
#    alphas = 2*mu.mat[j, ]
#    alphas = sapply(alphas, function(a) min(a, 1))
    points(rep(klist[j], nsim) + xjitter, pow.mat[j, ] + yjitter, pch = ".", col = rgb(red=0, green=0, blue=1), cex = 2)
  }
  dev.off()
}





