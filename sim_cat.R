
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_sim.R')
source('tex_table.R')


design = 'categorical'

nsim = 500
n = 100
num.nonzero = 10
k = num.nonzero
max.steps = 14
upper.coeff = 1.9
lower.coeff = 1.7

sigma = 1
groups = sort(rep(1:40, 5))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(g))
upper = upper.coeff*mult
lower = lower.coeff*mult

corr = 0
if ((corr != 0) & (design != 'gaussian')) {
  stop("nonzero only supported for gaussian design")
}

beta = beta_staircase(groups, num.nonzero, upper, lower, perturb = TRUE, cat.vars = 1:max(groups))
filename = paste0('figs/', design, '_size5_n', n, '_p', p, '_g', g, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
filename = paste0(filename, '.pdf')
pdf(filename)
output.l <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, design = design, corr = corr, categorical = TRUE, plot = TRUE)
dev.off()


groups = sort(rep(1:200, 5))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(g))
upper = upper.coeff*mult
lower = lower.coeff*mult
beta = beta_staircase(groups, num.nonzero, upper, lower, perturb = TRUE, cat.vars = 1:max(groups))
filename = paste0('figs/', design, '_size5_n', n, '_p', p, '_g', g, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
filename = paste0(filename, '.pdf')
pdf(filename)
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, design = design, corr = corr, categorical = TRUE, plot = TRUE)
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
