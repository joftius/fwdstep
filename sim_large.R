
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_sim.R')
source('tex_table.R')


design = 'ternary'
corr = 0 # nonzero only supported for gaussian design
nsim = 400
n = 100
num.nonzero = 8
max.steps = 12
mult = sqrt(2*log(p))
sigma = 1

groups = 1:200
p = length(groups)
upper.coeff = 1.8
lower.coeff = 1.7
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
filename = paste0(filename, '.pdf')
pdf(filename)
output.l <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, design = design, corr = corr, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(rep(1:30, 5), rep(31:35, 10)))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(g))
upper = upper.coeff*mult
lower = lower.coeff*mult
beta = beta_staircase(groups, num.nonzero, upper, lower)
filename = paste0('figs/', design, '_size5-10_n', n, '_p', p, '_g', g, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
filename = paste0(filename, '.pdf')
pdf(filename)
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, design = design, corr = corr, rand.beta = TRUE, plot = TRUE)
dev.off()


caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, m1))
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))

rownames(results.l) = paste("(1)", rownames(results.l))
rownames(results.g) = paste("(2)", rownames(results.g))

file = paste0(design, "_large_selection.tex")
tex_table(file, rbind(results.l, results.g), caption = caption)
