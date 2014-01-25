
# Use this as a template, it works in the fwd_step folder

source('../generate_data.R')
#source('fstep_minimal.R')


##########################
### Modify these lines ###
fname = "p20"           ##
ns = c(500) #, 200, 400)   ##
ps = 20 # floor(ns/2)        ##
nsim = 100               ##
klist = 1:6
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
  
  n = ns[i]
  p = ps[i]
  mult = sqrt(2*log(p+p*(p-1)/2))
  upper.coeff = 20.4
  lower.coeff = 10.4
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
      Y = X %*% beta + rnorm(n)
      Y = Y - mean(Y)

      # Weights?
      #weights = sqrt(rle(all.groups)$lengths)
      #weights = rep(1, max(all.groups)) #c(rep(1, p), rep(sqrt(2), max(all.groups) - p))
      
      # Fit forward stepwise
      added.groups = fstep_fit(X=X, Y=Y, groups=all.groups, max.steps=k)

      active.set = c()
      count = 0
      for (g in added.groups) {
        if (count < (4*k/3)) {
          active.set = c(active.set, g)
          count = count + ifelse(g <= p, 1, 3)
        }
      }
#      true.added = intersect(added.groups, true.groups)
#      false.added = setdiff(added.groups, true.groups)
      pow.mat[j, iter] = power_glinternet(k, true.ints, all.groups, int.groups, active.set)

    }
  }
  
  power.MCavg = rowMeans(pow.mat)

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





