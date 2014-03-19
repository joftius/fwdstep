
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_sim.R')
source('tex_table.R')
source('plots.R')

design = 'gaussian'
fn.append = ''
corr = 0 # nonzero only supported for gaussian design
noisecorr = 0

nsim = 100
n = 100
num.nonzero = 2
k = num.nonzero
max.steps = 8
upper.coeff = 0.8
lower.coeff = 0.5
Sigma = (1-noisecorr)*diag(rep(1,n)) + noisecorr
groups = 1:10
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(p))
upper = upper.coeff*mult
lower = lower.coeff*mult

if ((corr != 0) & (design != 'gaussian')) {
  stop("nonzero only supported for gaussian design")
}

beta = beta_staircase(groups, num.nonzero, upper, lower)

if (corr != 0) {
  fn.append = paste0(fn.append, '_corr', corr)
}
if (noisecorr != 0) {
  fn.append = paste0(fn.append, '_noisecorr', noisecorr)
}

print(c(upper, lower))

output.l <- fwd_group_simulation(n, Sigma, groups, beta, nsim, max.steps, design = design, corr = corr, rand.beta = TRUE)

with(output.l, step_plot(TrueStep, null.p, signal.p, chi.p, num.nonzero, n, p, g, ugsizes, max.steps, upper.coeff, lower.coeff, max.beta, min.beta, fwd.power, design, fn.append))


groups = sort(c(rep(1:8, 5), 9:10))
num.nonzero = 2
max.steps = 8
p = length(groups)
g = length(unique(groups))
mult = sqrt(2*log(p)) # or g?
upper = upper.coeff*mult
lower = lower.coeff*mult

beta = beta_staircase(groups, num.nonzero, upper, lower)

print(beta)
print(c(upper, lower))

output.g <- fwd_group_simulation(n, Sigma, groups, beta, nsim, max.steps, design = design, corr = corr, rand.beta = TRUE)

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
file = paste0(file, ".tex")

tex_table(file, rbind(results.l, results.g), caption = caption)
