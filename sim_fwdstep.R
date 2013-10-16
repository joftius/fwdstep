
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')
source('tex_table.R')

nsim = 20
n = 100
sigma = 1

groups = 1:200
p = length(groups)
mult = 1*sqrt(2*log(p))
upper = 1.5*mult
lower = 1.1*mult
num.nonzero = 8
max.steps = 12
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_lasso.pdf')
output.l <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(rep(1:20, 5), rep(21:30, 10)))
p = length(groups)
mult = 1*sqrt(2*log(p))
upper = 1.5*mult
lower = 1.1*mult
num.nonzero = 8
max.steps = 12
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_design.pdf')
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()


caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, m1))
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))

rownames(results.l) = paste("(1)", rownames(results.l))
rownames(results.g) = paste("(2)", rownames(results.g))

file = "large_sim_selection.tex"
tex_table(file, rbind(results.l, results.g), caption = caption)
