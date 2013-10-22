
source('fwd_step.R')
source('generate_data.R')
source('fwd_step_sim.R')
source('tex_table.R')

nsim = 25
sigma = 1
n = 100
p.list = c(25, 50, 200, 400)
k.list = c(3, 4, 5, 6, 7, 8, 9, 10)
max.steps = max(k.list) + 8
mu.mat = matrix(Inf, nrow = length(k.list), ncol = length(p.list))
mu.mat.sd = matrix(Inf, nrow = length(k.list), ncol = length(p.list))
pow.mat = matrix(0, nrow = length(k.list), ncol = length(p.list))
pow.mat.sd = matrix(0, nrow = length(k.list), ncol = length(p.list))

for (p in p.list) {
  j = which(p.list == p)
  groups = 1:p
  mult = 1*sqrt(2*log(p))
  upper = 1.2*mult
  lower = 1.0*mult  
  for (k in k.list) {
    i = which(k.list == k)
    beta = beta_staircase(groups, num.nonzero = k, upper, lower)
    output <- fwd_group_simulation(n, sigma, groups, beta, nsim, max.steps, rand.beta = TRUE, coherence = TRUE, plot = FALSE)
    mus = output$coherence
    pows = rowSums(output$true.step[, 1:k])/k
    mu.mat[i, j] = mean(mus)
    mu.mat.sd[i, j] = sd(mus)
    pow.mat[i, j] = mean(pows)
    pow.mat.sd[i, j] = sd(pows)
    cat(paste("Finished k =", k, "\n"))
  }
  cat(paste("Finished p =", p, "\n"))
}


#caption = "Evaluation of model selection using several stopping rules based on our p-values. The naive stopping rule performs well."
#results = with(output, sim_select_stats(signal.p, active.set, true.step, m1))
#file = "large_sim_selection____________what.tex"
#tex_table(file, results, caption = caption)
