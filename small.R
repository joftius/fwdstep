set.seed(1)
source("generate_data.R")
source("fwd_step.R")
source("tex_table.R")

n = 40
groups = c(1,1,1,2,2,3,3,3,3,4,4,5,5,6,6,7,7,7,8,8,9,9,10,10,10,10)
group.labels = unique(groups)
max.steps = 8
weights = rep(1, length(group.labels))
X = categorical_design(n, groups)
beta = rep(0, length(groups))
beta[groups == 1] = c(1,.5,-1)
beta[groups == 9] = c(.5,-.5)

Y = X %*% beta + rnorm(n)
Y = Y - mean(Y)
Sigma = diag(rep(1,n))
Xnormed = frob_normalize(X, groups)

results = forward_group(Xnormed, Y, groups, weights, Sigma, max.steps = max.steps, cat.groups = group.labels)

step = 1:max.steps
Variable = results$active.set
pval = results$p.vals
chi.pval = results$chi.pvals
matrix = rbind(Variable, pval, chi.pval)
colnames(matrix) = step
texfile = "tables/small.tex"
caption = "A small illustrative example"
tex_table(texfile, matrix, caption)
