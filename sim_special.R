#type = 'glint'
type = 'gamsel'


source('simulation.R')
source('tex_table.R')
source('plots.R')

estimation = FALSE

if (type == 'gamsel') {
    design = 'uniform'
    n = 500
    groups = 1:50
    nsim = 100
} else {
    design = 'gaussian'
    n = 100
    groups = 1:30
    nsim = 100
}

k = 6
k0 = 2
upper = 15
lower = 10
max.steps = 18


corr = 0 # nonzero only supported for gaussian design
noisecorr = 0

p = length(groups)

if ((corr != 0) & (design != 'gaussian')) {
    stop("Feature correlation only supported for gaussian design")
}

output = run_simulation(
    nsim = nsim,
    type = type,
    design = design,
    n = n, groups = groups, k = k, k0 = k0,
    upper = upper, lower = lower,
    max.steps = max.steps,
    estimation = estimation, verbose = TRUE)

special.power = paste0("(", round(output$special.fwd.power, 2), ")")
with(output, step_plot(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper, lower, max.beta, min.beta, fwd.power, design, filename, main.append = special.power))

ps.fname = paste0("figs/bysignal/", output$filename, ".pdf")
k = output$k
#psr = t(apply(output$psr.mat, 1, cumsum)/k)
psr = output$psr.mat
l = nrow(psr)
m = ncol(psr)
for (j in 1:nrow(psr)) {
    inds = which(psr[j,] >= 1)
      if (length(inds) > 1) {
            inds = inds[-1]
                psr[j,inds] = psr[j,inds] + cumsum(rep(1/(3*k), length(inds)))
          }
  }
psr = psr + matrix(0.02*rnorm(l*m), nrow=l)
pvals = output$signal.p
plot.main = paste0("n = ", n, ", p = ", p, ", signal strength ", lower, "/", upper)
pdf(ps.fname)
plot(psr[1, ], pvals[1, ], pch = ".", cex = 2, xlim=c(-0.1, max(psr) + .1), ylim=c(-0.1, 1.1), xlab = "Proportion of signal recovered", ylab = "P-values", main = plot.main)
for (j in 2:l) {
    inds = which(psr[j,] <= 2)
      points(psr[j, inds], pvals[j, inds], pch = ".", cex = 2)
  }
abline(h = c(.9, .7, .5, .3, .1), lty = 2, col = "gray")
abline(h = c(1, .8, .6, .4, .2, 0), lty = 3, col = "gray")
abline(v = 1, col = "gray")
dev.off()

#stop("No grouped features for GAMSel")

if (type == 'glint') {
  groups = sort(c(rep(1:15, 2), rep(16:20, 3)))
} else {
  n = 100
  groups = 1:100
  k = 3
  k0 = 1
}

p = length(groups)

output.g = run_simulation(
    nsim = nsim,
    type = type,
    design = design,
    n = n, groups = groups, k = k,
    upper = upper, lower = lower,
    max.steps = max.steps,
    estimation = estimation, verbose = TRUE)


with(output.g, step_plot(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper, lower, max.beta, min.beta, fwd.power, design, filename))


caption = "Evaluation of model selection using several stopping rules based on our p-values."
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, k))
results.l = with(output, sim_select_stats(signal.p, active.set, true.step, k))

file = paste0("tables/", output$filename, ".tex")

rownames(results.g) = paste("(g)", rownames(results.g))
tex_table(file, rbind(results.l, results.g), caption = caption)

if (estimation) {
    file = paste0("tables/", output$filename, "_estimation.tex")
    est.err = output$estimation
    est.err.g = output.g$estimation
    rownames(est.err.g) = paste("(g)", rownames(est.err.g))
    tex_table(file, rbind(est.err, est.err.g),
              caption = "Prediction and estimation errors")
}
