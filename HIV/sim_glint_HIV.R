#db="PI"
db="NRTI"

binary.encoding = TRUE
nsim = 100

setwd('..')
source('simulation.R')
source('generate_data.R')
source('glint/generate_glint.R')
source('tex_table.R')
source('plots.R')


if (db == "NRTI") {
    data = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/NRTI_DATA.txt", sep='\t', header=TRUE)
    resps = 4:9
} else {
    data = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt", sep = "\t", header = TRUE)
    resps = 4:10
}

type = "glint"
design = "fixed"

corr = 0 # nonzero only supported for gaussian design
noisecorr = 0
num.nonzero = 6
num.main = 2
k = num.nonzero
max.steps = 30
upper = 10
lower = 5

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
fixed.data = generate_glint(X=fixed.X, groups, cat.groups = groups)

## fixed.X = fixed.data$X
## all.groups = fixed.data$all.groups
## special.groups = fixed.data$special.groups
## default.groups = fixed.data$default.groups
## cat.groups = unique(groups)

n = nrow(fixed.X)
Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
p = length(groups)
g = length(unique(groups))

output = run_simulation(
    nsim = nsim,
    type = type,
    design = design,
    n = n, groups = groups, k = k,
    k0 = num.main,
    upper = upper, lower = lower,
    max.steps = max.steps,
    fixed.data = fixed.data,
    cat.groups = unique(groups),
    estimation = FALSE, verbose = TRUE
    )

output$filename = paste0("HIV_", output$filename)

special.power = paste0("(", round(output$special.fwd.power, 2), ")")

with(output, step_plot(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper, lower, max.beta, min.beta, fwd.power, design, filename, main.append = special.power))


ps.fname = paste0('figs/bysignal/', output$filename, ".pdf")
k = output.l$k
#psr = t(apply(output.l$psr.mat, 1, cumsum)/k)
psr = output.l$psr.mat
print(psr)
l = nrow(psr)
m = ncol(psr)
for (j in 1:nrow(psr)) {
    inds = which(psr[j,] >= 1)
      if (length(inds) > 1) {
            inds = inds[-1]
                psr[j,inds] = psr[j,inds] + cumsum(rep(1/(3*k), length(inds)))
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


caption = "Evaluation of model selection using several stopping rules based on our p-values."
results = with(output, sim_select_stats(signal.p, active.set, true.step, k))

file = paste0("tables/", output$filename, ".tex")

tex_table(file, results, caption = caption)
