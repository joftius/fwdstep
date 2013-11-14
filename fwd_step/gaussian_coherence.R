
source('grid_gaussian.R')
source('coherence.R')
source('../generate_data.R')
source('../tex_table.R')

filename = "coherence_gaussian.tex"
caption = "Mean coherence for Gaussian design matrices"

ns = c(100, 200, 400)
pmults = c(.5, 1, 2, 3)
nsim = 20

mu.means = matrix(NA, nrow=length(ns), ncol=length(pmults))
rownames(mu.means) = ns
colnames(mu.means) = pmults
mu.sds = mu.means

for (n in ns) {
  i = which(n == ns)
  
  for (pmult in pmults) {
    j = which(pmults == pmult)
    
    p = ceiling(n*pmult)
    groups = 1:p
    mus = c()

    if ((n > 200) & (pmult > 1)) {
      nsim = 5
    }

    for (iter in 1:nsim) {
      mat = gaussian_design(n, groups)
      mus = c(mus, coherence(mat))
    }

    mu.means[i, j] = mean(mus)
    mu.sds[i, j] = sd(mus)

    cat(as.character(Sys.time()), " --- Finished n=", n, "and p=", p, "\n")
  }
}

#f = file(filename)
mean.table = xtable(mu.means, caption = caption)
sd.table = xtable(mu.sds, caption = "Monte Carlo standard errors for above estimates")
cat(print(mean.table), file = filename)
cat(print(sd.table), file = filename, append=TRUE)
#close(f)



