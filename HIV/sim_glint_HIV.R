setwd('..')
source('glint/fwd_step.R')
source('generate_data.R')
source('glint/fwd_glint_sim.R')
#source('tex_table.R')
source('plots.R')

## PI = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/NRTI_DATA.txt", sep='\t', header=TRUE)
## db="PI"
## data=PI
## resps = 4:9
NRTI = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt", sep = "\t", header = TRUE)
data=NRTI
db="NRTI"
resps = 4:10

design = paste0("HIV_", db)
fn.append = ''
corr = 0 # nonzero only supported for gaussian design
noisecorr = 0
nsim = 20
num.nonzero = 3
k = num.nonzero
max.steps = 6
upper.coeff = 1.8
lower.coeff = 1.2

# Clean data, restrict to cleaned subset
nresps = length(resps)
Z = data[,resps]

clean_col = function(column) {
  outcol = as.character(column)
  outcol[column == '.'] = '-'
  return(as.factor(outcol))
}

cnames = c()
for (j in max(resps):ncol(data)) {
      newcol = clean_col(data[,j])
          if ((length(levels(newcol)) > 1) & (min(table(newcol)) >= 2)) {
                    cnames = c(cnames, names(data)[j])
                            Z = cbind(Z, newcol)
                  }
    }
colnames(Z) = c(names(data)[resps], cnames)
# Z should not be modified after this point

fixed.X = matrix(0, nrow=nrow(Z))
groups = c()
g = 0
for (j in (nresps+1):ncol(Z)) {
  g = g + 1
  submat = model.matrix(~ Z[,j] - 1)
  colnames(submat) = paste0(names(Z)[j], levels(Z[,j]))
  groups = c(groups, rep(g, ncol(submat)))
  fixed.X = cbind(fixed.X, submat)
}
fixed.X = fixed.X[,-1]
fixed.data = generate_glinternet(fixed.X, groups, cat.groups = groups)
fixed.X = fixed.data$X
groups = fixed.data$all.groups
int.groups = fixed.data$int.groups
cat.groups = unique(groups)
n = nrow(fixed.X)
Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
p = length(groups)
g = length(unique(groups))
weights = rep(1, g)
mult = sqrt(2*log(g))
upper = upper.coeff*mult
lower = lower.coeff*mult

if ((corr != 0) & (design != 'gaussian')) {
  stop("nonzero only supported for gaussian design")
}

#beta = beta_staircase(groups, num.nonzero, upper, lower)

if (corr != 0) {
  fn.append = paste0(fn.append, '_corr', corr)
}
if (noisecorr != 0) {
  fn.append = paste0(fn.append, '_noisecorr', noisecorr)
}
fn.append = paste0(fn.append, '_glint')

output.l <- fwd_glint_simulation(n, Sigma, groups, num.nonzero, lower, upper, nsim, max.steps, design = design, corr = corr, fixed.data = fixed.data, cat.groups = cat.groups)

with(output.l, step_plot(TrueStep, null.p, signal.p, chi.p, num.nonzero, n, p, g, ugsizes, max.steps, upper.coeff, lower.coeff, max.beta, min.beta, fwd.power, design, fn.append))

ps.fname = paste0('figs/bysignal/', design, '_size1_n', n, '_p', p, '_g', p, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  ps.fname = paste0(ps.fname, '_corr', corr)
}
if (noisecorr != 0) {
  ps.fname = paste0(ps.fname, '_noisecorr', noisecorr)
}
ps.fname = paste0(ps.fname, ".pdf")
m1 = output.l$m1
#psr = t(apply(output.l$psr.mat, 1, cumsum)/m1)
psr = output.l$psr.mat
print(psr)
l = nrow(psr)
m = ncol(psr)
for (j in 1:nrow(psr)) {
    inds = which(psr[j,] >= 1)
      if (length(inds) > 1) {
            inds = inds[-1]
                psr[j,inds] = psr[j,inds] + cumsum(rep(1/(3*m1), length(inds)))
          }
  }
psr = psr + matrix(0.04*rnorm(l*m), nrow=l)
pvals = output.l$signal.p
plot.main = paste0("n = ", n, ", p = ", p, ", signal strength ", lower.coeff, "/", upper.coeff)
pdf(ps.fname)
plot(psr[1, ], pvals[1, ], pch = ".", cex = 2, xlim=c(-0.1, max(psr) + .1), ylim=c(-0.1, 1.1), xlab = "Proportion of signal recovered", ylab = "P-values", main = plot.main)
for (j in 2:l) {
    inds = which(psr[j,] <= 2)
      points(psr[j, inds], pvals[j, inds], pch = ".", cex = 2)
  }
abline(h = c(.9, .7, .5, .3, .1), lty = 2, col = "gray")
abline(h = c(1, .8, .6, .4, .2, 0), lty = 3, col = "gray")
abline(v = 1, col = "gray")
dev.off()

stop("no selection")

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, m1))

#rownames(results.l) = paste("(1)", rownames(results.l))
rownames(results.g) = paste("(g)", rownames(results.g))

file = paste0("tables/", design, "_n", n, "_p", p, "_k", k, "_lower", lower.coeff, "_upper", upper.coeff)
if (corr != 0) {
  file = paste0(file, "_corr", corr)
}
file = paste0(file, ".tex")

tex_table(file, rbind(results.l, results.g), caption = caption)
