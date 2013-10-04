
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')
source('tex_table.R')

nsim = 500
n = 150
sigma = 1
groups = sort(c(rep(1:100, 5), rep(101:150, 10), 151:200))
p = length(groups)
upper = 3
lower = 2
num.nonzero = 10
max.steps = 14
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_design.pdf')
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()


groups = sort(c(groups, 151:200))
p = length(groups)
upper = 3 #(1+.1)*sqrt(2*log(p))
lower = 2 #(1-.1)*sqrt(2*log(p))
beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign = TRUE, perturb = TRUE, cat.vars = 1:max(groups))
warnings()

pdf('figs/large_sim_categorical_design.pdf')
output.c = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, categorical = TRUE, plot = TRUE)
dev.off()
warnings()

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))
results.c = with(output.c, sim_select_stats(signal.p, active.set, true.step, m1))
rownames(results.g) = paste("(1)", rownames(results.g))
rownames(results.c) = paste("(2)", rownames(results.c))
file = "large_sim_results.tex"
tex_table(file, rbind(results.g, results.c), caption = caption)
