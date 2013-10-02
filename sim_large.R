
library(xtable)

source('fwd_step.R')
source('generate_data.R')
source('fwd_step_test.R')

nsim = 100
n = 200
sigma = 1
groups = sort(c(rep(1:10, 5), rep(11:20, 10), rep(21:25, 15), 26:50))
p = length(groups)
upper = (1+.1)*sqrt(2*log(p))
lower = (1-.3)*sqrt(2*log(p))
num.nonzero = 10
max.steps = 14
beta = beta_staircase(groups, num.nonzero, upper, lower)

pdf('figs/large_sim_gaussian_design.pdf')
output.g <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, plot = TRUE)
dev.off()

#groups = sort(c(groups, 2, 8))
#beta = beta_staircase(groups, num.nonzero, upper, lower)
#pdf('figs/small_sim_categorical_design.pdf')
#output.c = fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, categorical = TRUE, plot = TRUE)
#dev.off()
#warnings()

caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
selection.results = with(output.g, xtable(sim_select_stats(signal.p, active.set, true.step, m1), caption = caption))
file = file("selection_sim_results.tex")
writeLines(print(selection.results), file)
close(file)


#with(output.c, print(sim_select_stats(signal.p, active.set, true.step, m1)))
