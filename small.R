set.seed(1)
source("generate_data.R")
source("fwd_step.R")
source("tex_table.R")

max.steps = 8
n = 40
groups = c(1,1,1,2,2,3,3,3,3,4,4,5,5,6,6,7,7,7,8,8,9,9,10,10,10,10)
group.labels = unique(groups)
cat.groups = group.labels
group.sizes = rle(groups)$lengths
weights = rep(1, length(group.labels))

X = categorical_design(n, groups)
beta = rep(0, length(groups))
beta[groups == 1] = c(1,.5,-1)
beta[groups == 9] = c(.5,-.5)

Y = X %*% beta + rnorm(n)
Y = Y - mean(Y)
Sigma = diag(rep(1,n))
Xnormed = frob_normalize(X, groups)

active.set = 0
eff.p = 0
p.vals = c()
chi.pvals = c()
MC.pvals = c()
c.vars = c()
Ls = c()
ranks = c()
Ynorms = c()

MC_pval = function(L, X, inactive.groups, inactive.labels, inactive.inds, M=200) {
  n = nrow(X)
  Ls = c()
  for (iter in 1:M) {
    epsilon = rnorm(n)
    U = t(X[, inactive.inds]) %*% epsilon
    terms = sapply(inactive.labels, function(x) sqrt(sum(U[inactive.groups==x]^2)))
    Ls = c(Ls, max(terms))
  }
  q = sum(Ls > L)/M
  return(q)
}

Y.update = Y
X.update = X
MC.time = 0
Tx.time = 0

for (i in 1:max.steps) {
  # Tchi statistic computation
  before.Tx.time = as.numeric(Sys.time())
  output = add_group(X.update, Y.update, groups, weights, Sigma, active.set, eff.p, cat.groups = cat.groups)
  Tx.time = Tx.time + as.numeric(Sys.time()) - before.Tx.time
  L = output$test.output[[1]]
  inactive.labels = setdiff(group.labels, active.set)
  inactive.inds = sapply(groups, function(x) is.element(x, inactive.labels))
  inactive.groups = groups[inactive.inds]
  # MC statistic computation
  before.MC.time = as.numeric(Sys.time())
  MC.pval = MC_pval(L, X, inactive.groups, inactive.labels, inactive.inds)
  MC.time = MC.time + as.numeric(Sys.time()) - before.MC.time
  active.set = output$active.set
  imax = output$added
  grank = output$grank
  ranks = c(ranks, grank)
  eff.p = output$eff.p
  RSS = sum(Y.update^2)
  Y.update = output$Y.update
  X.update = output$X.update
  Ynorms = c(Ynorms, RSS)    
  RSSdrop = RSS - sum(Y.update^2)
  p.vals = c(p.vals, output$p.value)
  chi.p = pchisq(RSSdrop, lower.tail=F, df=grank)
  chi.pvals = c(chi.pvals, chi.p)
  MC.pvals = c(MC.pvals, MC.pval)
  c.vars = c(c.vars, output$var)
  Ls = c(Ls, output$test.output[[1]])
}

step = 1:max.steps
Variable = active.set
pval = p.vals
chi.pval = chi.pvals
MC.pval = MC.pvals
matrix = rbind(Variable, pval, chi.pval, MC.pval)
rownames(matrix) = c("Variable", "Tchi", "chi2", "max-chi")
colnames(matrix) = step
texfile = "tables/small.tex"
caption = "A small illustrative example. Elapsed time for forward stepwise textitand $T\\chi$ computation:"
caption = paste(caption, round(Tx.time, 3))
caption = paste(caption, "seconds. Elapsed time for 200 Monte Carlo sample estimate of max-chi:")
caption = paste(caption, round(MC.time, 3))
caption = paste(caption, "seconds.")
tex_table(texfile, matrix, caption)
