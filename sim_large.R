
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')
source('tex_table.R')

nsim = 200
n = 150
sigma = 1

groups = 1:1050
p = length(groups)
mult = sqrt(2*log(p))
upper = 1.1*mult
lower = 0.5*mult
num.nonzero = 20
max.steps = 30
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_lasso.pdf')
output.l <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(rep(1:100, 5), rep(101:150, 10), 151:200))
p = length(groups)
mult = sqrt(2*log(p))
upper = 1.5*mult
lower = 1.1*mult
num.nonzero = 8
max.steps = 12
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_design.pdf')
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(groups, 151:200))
p = length(groups)
mult = sqrt(2*log(p))
upper = 1.5*mult
lower = 1.1*mult
beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign = TRUE, perturb = TRUE, cat.vars = 1:max(groups))
warnings()

pdf('figs/large_sim_categorical_design.pdf')
output.c = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, categorical = TRUE, plot = TRUE)
dev.off()
warnings()

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.l = with(output.l, sim_select_stats(signal.p, active.set, true.step, m1))
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))
results.c = with(output.c, sim_select_stats(signal.p, active.set, true.step, m1))
rownames(results.l) = paste("(1)", rownames(results.l))
rownames(results.g) = paste("(2)", rownames(results.g))
rownames(results.c) = paste("(3)", rownames(results.c))
file = "large_sim_results.tex"
tex_table(file, rbind(results.g, results.c), caption = caption)
