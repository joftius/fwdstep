
source('fwd_step.R')
source('../generate_data.R')
source('fwd_glint_sim.R')
#source('tex_table.R')


design = 'gaussian'
corr = 0 # nonzero only supported for gaussian design

nsim = 200
n = 250
num.nonzero = 6
k = num.nonzero
max.steps = 10
upper.coeff = 8
lower.coeff = 5

sigma = 1
groups = 1:20
p = length(groups)
mult = sqrt(2*log(p))
upper = upper.coeff*mult
lower = lower.coeff*mult

if ((corr != 0) & (design != 'gaussian')) {
  stop("nonzero only supported for gaussian design")
}

beta = beta_staircase(groups, num.nonzero, upper, lower)
filename = paste0('../figs/', design, '_size1_n', n, '_p', p, '_g', p, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
filename = paste0(filename, '_glint.pdf')

pdf(filename)
output.l <- fwd_glint_simulation(n, sigma, groups, num.nonzero, lower, upper, nsim, max.steps, design = design, corr = corr, plot = TRUE)
dev.off()

ps.fname = paste0('../figs/bysignal/', design, '_size1_n', n, '_p', p, '_g', p, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  ps.fname = paste0(ps.fname, '_corr', corr)
}
ps.fname = paste0(ps.fname, '_glint.pdf')
m1 = output.l$m1
#psr = t(apply(output.l$psr.mat, 1, cumsum)/m1)
psr = output.l$psr.mat
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



#stop("Not this time!")

num.nonzero = 6
k = 6
max.steps = 9
groups = sort(c(rep(1:15, 2), rep(16:20, 3)))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(g))
upper = upper.coeff*mult
lower = lower.coeff*mult

filename = paste0('../figs/', design, '_size2-3_n', n, '_p', p, '_g', g, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
filename = paste0(filename, '_glint.pdf')
pdf(filename)
output.g <- fwd_glint_simulation(n, sigma, groups, num.nonzero, lower, upper, nsim, max.steps, design = design, corr = corr, plot = TRUE)
dev.off()

print(warnings())

print(summary(output.l$main1))
print(summary(output.l$int37))

print(summary(output.g$main1))
print(summary(output.g$int37))


stop("No selection stats for this")

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, m1))
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))

#rownames(results.l) = paste("(1)", rownames(results.l))
rownames(results.g) = paste("(g)", rownames(results.g))

file = paste0("tables/", design, "_n", n, "_p", p, "_k", k, "_lower", lower.coeff, "_upper", upper.coeff)
if (corr != 0) {
  file = paste0(file, "_corr", corr)
}
file = paste0(file, ".tex")

tex_table(file, rbind(results.l, results.g), caption = caption)
