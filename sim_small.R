
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_sim.R')
source('tex_table.R')


design = 'gaussian'
corr = 0 # nonzero only supported for gaussian design
noisecorr = .1

nsim = 200
n = 100
num.nonzero = 5
k = num.nonzero
max.steps = 10
upper.coeff = 1.2
lower.coeff = 1.1
Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
groups = 1:50
p = length(groups)
mult = sqrt(2*log(p))
upper = upper.coeff*mult
lower = lower.coeff*mult

if ((corr != 0) & (design != 'gaussian')) {
  stop("nonzero only supported for gaussian design")
}

beta = beta_staircase(groups, num.nonzero, upper, lower)
filename = paste0('figs/', design, '_size1_n', n, '_p', p, '_g', p, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
if (noisecorr != 0) {
  filename = paste0(filename, '_noisecorr', noisecorr)
}
filename = paste0(filename, '.pdf')
print(c(upper, lower))
pdf(filename)
output.l <- fwd_group_simulation(n, Sigma, groups, beta, nsim, max.steps, design = design, corr = corr, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(rep(1:8, 5), 9:18))
num.nonzero = 3
max.steps = 9
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(p)) # or g?
upper = upper.coeff*mult
lower = lower.coeff*mult

beta = beta_staircase(groups, num.nonzero, upper, lower)
filename = paste0('figs/', design, '_size1-5_n', n, '_p', p, '_g', g, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
if (noisecorr != 0) {
  filename = paste0(filename, '_noisecorr', noisecorr)
}
filename = paste0(filename, '.pdf')
print(beta)
print(c(upper, lower))
pdf(filename)
output.g <- fwd_group_simulation(n, Sigma, groups, beta, nsim, max.steps, design = design, corr = corr, rand.beta = TRUE, plot = TRUE)
dev.off()


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
