
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')
source('tex_table.R')

nsim = 100
n = 100
sigma = 1
groups = c(1, 1, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9, 9, rep(10, 10))
p = length(groups)
mult = 1*sqrt(2*log(p))
upper = 1.5*mult
lower = 1.1*mult
num.nonzero = 3
max.steps = 8
beta = beta_staircase(groups, num.nonzero, upper, lower)
beta_old = beta

pdf('figs/small_sim_gaussian_design.pdf')
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, predictions = TRUE, plot = TRUE)
dev.off()
warnings()


groups = sort(c(groups, 2, 8))
p = length(groups)
mult = 1*sqrt(2*log(p))
upper = 1.5*mult
lower = 1.1*mult
beta = beta_staircase(groups, num.nonzero, upper, lower, rand.sign = TRUE, perturb = TRUE, cat.vars = 1:max(groups))
warnings()

pdf('figs/small_sim_categorical_design.pdf')
output.c = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, categorical = TRUE, plot = TRUE)
dev.off()
warnings()

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, m1))
results.c = with(output.c, sim_select_stats(signal.p, active.set, true.step, m1))
rownames(results.g) = paste("(1)", rownames(results.g))
rownames(results.c) = paste("(2)", rownames(results.c))
file = "small_sim_selection.tex"
tex_table(file, rbind(results.g, results.c), caption = caption)


caption = "Prediction and estimation errors for a small simulation"
p.results.g = output.g$pred.err
file = "small_sim_estimation.tex"
tex_table(file, p.results.g, caption = caption)


#print(c(max(abs(beta)), max(abs(beta_old)), upper, lower))
#print(beta_old)
#print(beta)
