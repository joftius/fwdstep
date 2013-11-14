
source('grid_gaussian.R')
source('coherence.R')
source('../generate_data.R')
source('../tex_table.R')

filename = "coherence_gaussian_corr.tex"
caption = "Mean coherence for Gaussian design matrices with equi-correlation"

n = 100
rs = c(.01, .05, .1, .2)
pmults = c(.5, 1, 2, 3)
nsim = 20

mu.means = matrix(NA, nrow=length(rs), ncol=length(pmults))
rownames(mu.means) = rs
colnames(mu.means) = pmults
mu.sds = mu.means

for (r in rs) {
  i = which(r == rs)
  
  for (pmult in pmults) {
    j = which(pmults == pmult)
    
    p = ceiling(n*pmult)
    groups = 1:p
    mus = c()

    for (iter in 1:nsim) {
      mat = gaussian_design(n, groups, corr = r)
      mus = c(mus, coherence(mat))
    }

    mu.means[i, j] = mean(mus)
    mu.sds[i, j] = sd(mus)

    cat(as.character(Sys.time()), " --- Finished r=", r, "and p=", p, "\n")
  }
}

mean.table = xtable(mu.means, caption = caption)
sd.table = xtable(mu.sds, caption = "Monte Carlo standard errors for above estimates")
cat(print(mean.table), file = filename)
cat(print(sd.table), file = filename, append=TRUE)

