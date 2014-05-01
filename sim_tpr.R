
source('tpr_simulation.R')

#groups = 1:50
groups = 1:200
#groups = sort(c(rep(1:26, 2), rep(27:38, 4)))
#groups = sort(rep(1:36, 2))
#groups = sort(rep(1:50, 3))
#groups = sort(rep(1:200, 2))

type = 'default'
design = 'gaussian'
#design = 'categorical'
#design = 'ternary'
#design = 'uniform'

corr = .1
noisecorr = 0
#corr = 0
#noisecorr = 0

nsim = 100
n = 100

kmax = 24
upper = 2
lower = 1.2

### No editing below ###

if (type == 'gamsel') {
  design = 'uniform'
}

if (design == 'categorical') {
  cat.groups = unique(groups)
} else {
  cat.groups = NULL
}

if (design != 'gaussian') {
  corr = 0
}

estimation = FALSE


if ((corr != 0) & (design != 'gaussian')) {
    stop("Feature correlation only supported for gaussian design")
}

output = run_tpr_simulation(
    nsim = nsim,
    type = type,
    design = design,
    n = n, groups = groups, kmax = kmax,
    upper = upper, lower = lower,
    cat.groups = cat.groups,
    noisecorr = noisecorr,
    corr = corr,
    estimation = estimation, verbose = TRUE)


