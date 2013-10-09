
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')
source('tex_table.R')

nsim = 200
n = 150
sigma = 1

groups = 1:500
p = length(groups)
mult = sqrt(2*log(p))
upper = 1.1*mult
lower = 0.9*mult
num.nonzero = 20
max.steps = 30
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_lasso.pdf')
output.l <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(rep(1:50, 5), rep(51:75, 10)))
p = length(groups)
mult = sqrt(2*log(p))
upper = 1.1*mult
lower = 0.9*mult
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
