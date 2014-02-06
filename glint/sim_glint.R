
source('fwd_step.R')
source('../generate_data.R')
source('fwd_glint_sim.R')
#source('tex_table.R')


design = 'gaussian'
corr = 0 # nonzero only supported for gaussian design

nsim = 150
n = 200
num.nonzero = 6
k = num.nonzero
max.steps = 12
upper.coeff = 1.9
lower.coeff = 1.1

sigma = 1
groups = 1:20
p = length(groups)
mult = sqrt(2*log(p*(p+1)/2))
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

stop("no")

num.nonzero = 6
k = 6
max.steps = 9
groups = sort(c(rep(1:15, 2), rep(16:20, 4)))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(g))
upper = upper.coeff*mult
lower = lower.coeff*mult

filename = paste0('../figs/', design, '_size2-4_n', n, '_p', p, '_g', g, '_k', num.nonzero, '_lower', lower.coeff, '_upper', upper.coeff)
if (corr != 0) {
  filename = paste0(filename, '_corr', corr)
}
filename = paste0(filename, '_glint.pdf')
pdf(filename)
output.g <- fwd_glint_simulation(n, sigma, groups, num.nonzero, lower, upper, nsim, max.steps, design = design, corr = corr, plot = TRUE)
dev.off()

print(warnings())

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
