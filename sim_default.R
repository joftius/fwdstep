source('simulation.R')
source('tpr_simulation.R')
source('tex_table.R')
source('plots.R')

nsim = 400
type = 'default'
design = 'gaussian'
estimation = FALSE #TRUE
n = 100
#groups = 1:500
#groups = 1:400
groups = 1:200
#groups = 1:100
#groups = 1:50

k = 10
upper = 3
lower = 2

max.steps = 15
#max.steps = min(n/2, max(groups)-1) #floor(max(groups)/3)
corr = 0.2 # nonzero only supported for gaussian design
noisecorr = 0.1

if ((corr != 0) & (design != 'gaussian')) {
    stop("Feature correlation only supported for gaussian design")
}

output = run_simulation(
    nsim = nsim,
    type = type,
    design = design,
    n = n, groups = groups, k = k,
    upper = upper, lower = lower,
    max.steps = max.steps,
    corr = corr,
    noisecorr = noisecorr,
    estimation = estimation, verbose = TRUE)

output$main.append = paste0(", RIC: (", round(output$ric.tpp, 2), "/", output$ric.size, ")")

with(output, step_plot(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper, lower, max.beta, min.beta, fwd.power, design, filename, main.append))

## ps.fname = paste0('figs/bysignal/', design, '_size1_n', n, '_p', p, '_g', p, '_k', k, '_lower', lower, '_upper', upper)
## if (corr != 0) {
##   ps.fname = paste0(ps.fname, '_corr', corr)
## }
## if (noisecorr != 0) {
##   ps.fname = paste0(ps.fname, '_noisecorr', noisecorr)
## }
## ps.fname = paste0(ps.fname, '.pdf')
## k = output$k
## #psr = t(apply(output$psr.mat, 1, cumsum)/k)
## psr = output$psr.mat
## print(psr)
## l = nrow(psr)
## m = ncol(psr)
## for (j in 1:nrow(psr)) {
##     inds = which(psr[j,] >= 1)
##       if (length(inds) > 1) {
##             inds = inds[-1]
##                 psr[j,inds] = psr[j,inds] + cumsum(rep(1/(3*k), length(inds)))
##           }
##   }
## psr = psr + matrix(0.02*rnorm(l*m), nrow=l)
## pvals = output$signal.p
## plot.main = paste0("n = ", n, ", p = ", p, ", signal strength ", lower, "/", upper)
## pdf(ps.fname)
## plot(psr[1, ], pvals[1, ], pch = ".", cex = 2, xlim=c(-0.1, max(psr) + .1), ylim=c(-0.1, 1.1), xlab = "Proportion of signal recovered", ylab = "P-values", main = plot.main)
## for (j in 2:l) {
##     inds = which(psr[j,] <= 2)
##       points(psr[j, inds], pvals[j, inds], pch = ".", cex = 2)
##   }
## abline(h = c(.9, .7, .5, .3, .1), lty = 2, col = "gray")
## abline(h = c(1, .8, .6, .4, .2, 0), lty = 3, col = "gray")
## abline(v = 1, col = "gray")
## dev.off()

groups = sort(c(rep(1:30, 5), rep(31:35, 10)))
#groups = sort(c(rep(1:10, 3), rep(11:20, 2)))
#groups = 1:200
corr = 0
noisecorr = 0
max.steps = 10
k = 5

output.g = run_simulation(
    nsim = nsim,
    type = type,
    design = design,
    n = n, groups = groups, k = k,
    upper = upper, lower = lower,
    max.steps = max.steps,
    corr = corr,
    noisecorr = noisecorr,
    estimation = estimation, verbose = TRUE)

#output.g$main.append = paste0(", RIC: (", round(output.g$ric.tpp, 2), "/", output.g$ric.size, ")")

with(output.g, step_plot(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper, lower, max.beta, min.beta, fwd.power, design, filename))

caption = "Evaluation of model selection using several stopping rules based on our p-values."
results.l = with(output, sim_select_stats(signal.p, active.set, true.step, k, bic.list = bic.list, ric.list = ric.list))
results.g = with(output.g, sim_select_stats(signal.p, active.set, true.step, k))

#rownames(results.l) = paste0("*", rownames(results.l))
#rownames(results.g) = paste0("**", rownames(results.g))

file = paste0("tables/", output$filename, ".tex")

tex_table(file, rbind(results.l, results.g), caption = caption)

if (estimation) {
    file = paste0("tables/", output$filename, "_estimation.tex")
    est.err = output$estimation
    est.err.g = output.g$estimation
    rownames(est.err) = paste("(5)", rownames(est.err))
    rownames(est.err.g) = paste("(10)", rownames(est.err.g))
    tex_table(file, rbind(est.err, est.err.g),
              caption = "Prediction and estimation errors")
}
