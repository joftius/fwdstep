db = "NRTI"
#db = "PI"

binary.encoding = FALSE

setwd('..')
source('fwd_step.R')
source('generate_data.R')
source('simulation.R')
source('tex_table.R')
source('plots.R')
library(plyr)

if (db == "NRTI") {
  data = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/NRTI_DATA.txt", sep='\t', header=TRUE)
  resps = 4:9
} else {
  data = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt", sep = "\t", header = TRUE)
  resps = 4:10
}

nsim = 400
type = "default"
design = "fixed"
fn.append = '' #'Pg'
corr = 0 # nonzero only supported for gaussian design
noisecorr = 0
num.nonzero = 5
k = num.nonzero
max.steps = 10
upper = 8
lower = 4

# Clean data, restrict to cleaned subset
nresps = length(resps)
Z = data[,resps]

clean_col = function(column) {
  outcol = as.character(column)
  outcol[column == '.'] = '-'
  return(as.factor(outcol))
}

cnames = c()
for (j in (max(resps)+1):ncol(data)) {
  newcol = clean_col(data[,j])
  if (length(unique(newcol)) > 1) {
    newcol = mapvalues(newcol, from=levels(newcol), to=1:length(levels(newcol)))
    cnames = c(cnames, names(data)[j])
    Z = cbind(Z, newcol)
  }
}
colnames(Z) = c(names(data)[resps], cnames)

### Z should not be modified after this point ###
  
fixed.X = matrix(NA, nrow=nrow(Z))
groups = c()
g = 0
for (j in (nresps+1):ncol(Z)) {
  g = g + 1
  submat = model.matrix(~ Z[,j] - 1)
  if (binary.encoding) {
    if (ncol(submat) > 2) {
      submat = cbind(submat[,1], rowSums(submat[,2:ncol(submat)]))
    }
    colnames(submat) = paste0(names(Z)[j], c("-", "*"))
  } else {
    colnames(submat) = paste0(names(Z)[j], levels(Z[,j]))
  }
  groups = c(groups, rep(g, ncol(submat)))
  fixed.X = cbind(fixed.X, submat)
}
fixed.X = fixed.X[,-1]
cat.groups = unique(groups)
X.cat = Z[, (nresps+1):ncol(Z)]
colnames(X.cat) = paste0("X", cat.groups)
fixed.data = list(X=fixed.X, X.cat=X.cat)

n = nrow(fixed.X)
Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
p = length(groups)

output.l = run_simulation(
  nsim = nsim,
  type = type, design = design,
  n = n, groups = groups, k = k,
  upper = upper, lower = lower,
  max.steps = max.steps,
  cat.groups = cat.groups,
  fixed.data = fixed.data,
  estimation = FALSE, verbose = TRUE)


with(output.l, step_plot(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper, lower, max.beta, min.beta, fwd.power, design, filename))

## ps.fname = paste0('figs/bysignal/', design, '_size1_n', n, '_p', p, '_g', p, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
## if (corr != 0) {
##   ps.fname = paste0(ps.fname, '_corr', corr)
## }
## if (noisecorr != 0) {
##   ps.fname = paste0(ps.fname, '_noisecorr', noisecorr)
## }
## ps.fname = paste0(ps.fname, ".pdf")
## m1 = output.l$m1
## #psr = t(apply(output.l$psr.mat, 1, cumsum)/m1)
## psr = output.l$psr.mat
## print(psr)
## l = nrow(psr)
## m = ncol(psr)
## for (j in 1:nrow(psr)) {
##     inds = which(psr[j,] >= 1)
##       if (length(inds) > 1) {
##             inds = inds[-1]
##                 psr[j,inds] = psr[j,inds] + cumsum(rep(1/(3*m1), length(inds)))
##           }
##   }
## psr = psr + matrix(0.04*rnorm(l*m), nrow=l)
## pvals = output.l$signal.p
## plot.main = paste0("n = ", n, ", p = ", p, ", signal strength ", lower.coeff, "/", upper.coeff)
## pdf(ps.fname)
## plot(psr[1, ], pvals[1, ], pch = ".", cex = 2, xlim=c(-0.1, max(psr) + .1), ylim=c(-0.1, 1.1), xlab = "Proportion of signal recovered", ylab = "P-values", main = plot.main)
## for (j in 2:l) {
##     inds = which(psr[j,] <= 2)
##       points(psr[j, inds], pvals[j, inds], pch = ".", cex = 2)
##   }
## abline(h = c(.9, .7, .5, .3, .1), lty = 2, col = "gray")
## abline(h = c(1, .8, .6, .4, .2, 0), lty = 3, col = "gray")
## abline(v = 1, col = "gray")
## dev.off()

caption = "Evaluation of model selection using several stopping rules based on our p-values."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, k))

#rownames(results.l) = paste("(1)", rownames(results.l))
#rownames(results.g) = paste("(g)", rownames(results.g))

file = paste0("tables/", output.l$filename, ".tex")

tex_table(file, results.l, caption = caption)
