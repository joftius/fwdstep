
source('../generate_data.R')
source('fstep_minimal.R')
source('grid_gaussian.R')
source('coherence.R')

##########################
### Modify these lines ###
#fname = "p.5n"           ##
#ns = c(200) #, 200, 400)   ##
#ps = floor(ns/2)        ##
#nsim = 60               ##
#kplus = 20              ##
#########################
fname = "p4n"          ##
ns = c(100) #, 200)    ##
ps = ceiling(ns*4)     ##
nsim = 50              ##
kplus = 8              ##
##########################
mus = n_to_mu(ns)
ks = mu_to_k(mus)

for (i in 1:length(ns)) {
  mu.mat = matrix(0, nrow=1+kplus, ncol = nsim)
  pow.mat = mu.mat
  
  n = ns[i]
  p = ps[i]
  mult = sqrt(2*log(p))
  upper.coeff = 1.4
  lower.coeff = 1.4
  upper = upper.coeff*mult
  lower = lower.coeff*mult
  mu = mus[i]
  kmin = 1 # ks[i]
  kmax = kmin + kplus
  # Something
  klist = 0:(kmax-1)
  klist = 2*klist + 1
  groups = 1:p

  for (j in 1:(1+kplus)) {
    k = klist[j]
    print(k)
    
    for (iter in 1:nsim) {
      # Generate data
      beta = beta_staircase(groups, k, upper, lower, rand.sign=TRUE, permute=TRUE)
      true.groups = true_active_groups(groups, beta)
      X = gaussian_design(n, groups)
      Y = X %*% beta + rnorm(n)
      
      # Calculate coherence
      mu.mat[j, iter] = coherence(X)
      
      # Fit forward stepwise
      added.groups = fstep_fit(X, Y, groups, k)
      pow.mat[j, iter] = length(intersect(added.groups, true.groups)) / k
    }
  }

  power.MCavg = rowMeans(pow.mat)

  print(pow.mat)

  print(mu.mat)

  # Save plot
  plot.main = paste0("n = ", n, ", p = ", p, ", signal strength ", lower.coeff, "/", upper.coeff)
  pdf(paste0('../figs/fwdstep/gaussian_', fname, "_n", n, "_lower", lower.coeff, "_upper", upper.coeff, ".pdf"))
  plot(klist, power.MCavg, type = "l", main = plot.main, xlab = "Sparsity", ylab = "Average power", ylim = c(-0.1, 1.1), lwd = 2)
  for (j in 1:(1+kplus)) {
    xjitter = 0.05*rnorm(nsim)
    yjitter = 0.01*rnorm(nsim)
    alphas = 2*mu.mat[j, ]
    alphas = sapply(alphas, function(a) min(a, 1))
    points(rep(klist[j], nsim) + xjitter, pow.mat[j, ] + yjitter, pch = ".", col = rgb(red=0, green=0, blue=1, alpha = alphas), cex = 2)
  }
  dev.off()
}





