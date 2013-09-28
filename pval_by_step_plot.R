
library(xtable)

source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')

nsim = 500
n = 50
sigma = 1
groups = c(1, 1, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9, 9, rep(10, 10))
p = length(groups)
upper = (1+.1)*sqrt(2*log(p))
lower = (1-.1)*sqrt(2*log(p))
num.nonzero = 3
max.steps = 8
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/small_sim_gaussian_design.pdf')
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, plot = TRUE)
dev.off()

groups = sort(c(groups, 2, 8))
beta = beta_staircase(groups, num.nonzero, upper, lower)

#pdf('figs/small_sim_categorical_design.pdf')
#output.c = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, categorical = TRUE, plot = TRUE)
#dev.off()
#warnings()

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
selection.results = with(output.g, xtable(sim_select_stats(signal.p, active.set, true.step, m1), caption = caption))
file = file("small_sim_results.tex")
writeLines(print(selection.results), file)
close(file)


#with(output.c, print(sim_select_stats(signal.p, active.set, true.step, m1)))
