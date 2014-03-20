
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_sim.R')
source('tex_table.R')
source('plots.R')

design = 'categorical'
fn.append = ''
corr = 0
noisecorr = 0

nsim = 400
n = 100
num.nonzero = 4
k = num.nonzero
max.steps = 10
upper.coeff = 5
lower.coeff = 4
Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
groups = sort(rep(1:20, 5))
#groups = sort(rep(1:20, 2))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(p)/n)
upper = upper.coeff*mult
lower = lower.coeff*mult

if ((corr != 0) & (design != 'gaussian')) {
  stop("nonzero only supported for gaussian design")
}
if (noisecorr != 0) {
  fn.append = paste0(fn.append, "_noisecorr", noisecorr)
}

beta = beta_staircase(groups, num.nonzero, upper, lower, permute = TRUE, perturb = TRUE, cat.groups = unique(groups))

output.l = fwd_group_simulation(n, Sigma, groups, beta, nsim, max.steps, design = design, corr = corr, categorical = TRUE, rand.beta = TRUE, cat.groups = unique(groups))

print(warnings())

with(output.l, step_plot(TrueStep, null.p, signal.p, chi.p, num.nonzero, n, p, g, ugsizes, max.steps, upper.coeff, lower.coeff, max.beta, min.beta, fwd.power, design, fn.append))

groups = sort(c(rep(1:100, 4), rep(101:200, 2)))
#groups = sort(rep(1:30, 5))
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(p)/n)
upper = upper.coeff*mult
lower = lower.coeff*mult
beta = beta_staircase(groups, num.nonzero, upper, lower, perturb = TRUE, cat.groups = 1:max(groups))

output.g = fwd_group_simulation(n, Sigma, groups, beta, nsim, max.steps, design = design, corr = corr, categorical = TRUE, rand.beta = TRUE, cat.groups = unique(groups))

with(output.g, step_plot(TrueStep, null.p, signal.p, chi.p, num.nonzero, n, p, g, ugsizes, max.steps, upper.coeff, lower.coeff, max.beta, min.beta, fwd.power, design, fn.append))


caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, m1))
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))

#rownames(results.l) = paste("(1)", rownames(results.l))
rownames(results.g) = paste("(g)", rownames(results.g))

file = paste0("tables/", design, "_n", n, "_p", p, "_k", k, "_lower", lower.coeff, "_upper", upper.coeff)
if (corr != 0) {
  file = paste0(file, "_corr", corr)
}
if (noisecorr != 0) {
  file = paste0(file, "_noisecorr", noisecorr)
}
file = paste0(file, ".tex")

tex_table(file, rbind(results.l, results.g), caption = caption)
