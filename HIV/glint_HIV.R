setwd('..')
source('glint/fwd_step.R')
source('generate_data.R')
source('stop_rules.R')
source('tex_table.R')

PI = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/NRTI_DATA.txt", sep='\t', header=TRUE)
db="PI"
data=PI
resps = 4:9
## NRTI = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt", sep = "\t", header = TRUE)
## data=NRTI
## db="NRTI"
## resps = 4:10

nresps = length(resps)
max.steps = 20

# Clean data, restrict to cleaned subset
Z = data[,resps]

clean_col = function(column) {
    outcol = as.character(column)
    outcol[column == '.'] = '-'
    return(as.factor(outcol))
}

cnames = c()
for (j in max(resps):ncol(data)) {
    newcol = clean_col(data[,j])
    if ((length(levels(newcol)) > 1) & (min(table(newcol)) >= 2)) {
        cnames = c(cnames, names(data)[j])
        Z = cbind(Z, newcol)
    }
}
colnames(Z) = c(names(data)[resps], cnames)
# Z should not be modified after this point

## fs.data = matrix(0, nrow=nrow(Z))
## groups = c()
## g = 1
## for (j in (nresps+1):ncol(Z)) {
##     submat = model.matrix(~ Z[,j] - 1)
##     colnames(submat) = paste0(names(Z)[j], levels(Z[,j]))
##     groups = c(groups, rep(g, ncol(submat)))
##     g = g + 1
##     fs.data = cbind(fs.data, submat)
## }
## fs.data = fs.data[,-1]

# Perform analyses for each possible response

A.sets = list()
ns = c()
Gs = c()
ps = c()
Pvals = matrix(NA, nrow=nresps, ncol=max.steps)
Pvals.chi = Pvals
for (respind in 1:nresps) {
    subrows = which(!is.na(Z[,respind]))
    Y = log(Z[subrows,respind])
    n = length(Y)
    # Now restrict to gene variants with >= 2 observations
    X = matrix(NA, nrow=n)
    groups = c()
    g = 1
    cnames = c()
    for (j in (nresps+1):ncol(Z)) {
        newcol = Z[subrows, j]
        if ((length(levels(newcol)) > 1) & (min(table(newcol)) >= 2)) {
            cnames = c(cnames, names(Z)[j])
            submat = model.matrix(~ Z[subrows,j] - 1)
            colnames(submat) = paste0(names(Z)[j], levels(Z[,j]))
            groups = c(groups, rep(g, ncol(submat)))
            g = g + 1
            X = cbind(X, submat)
        }
    }
    X = X[ , -1]
    gld = generate_glinternet(X, groups)
    glX = gld$X
    all.groups = gld$all.groups
    glX = frob_normalize(glX, all.groups)
    G = max(all.groups)
    Gs = c(Gs, G)
    ns = c(ns, n)
    ps = c(ps, length(all.groups))
    #weights = sqrt(rle(groups)$lengths)
    weights = rep(1, max(all.groups))
    Sigma = diag(rep(1,n))
    # Perform forward stepwise
    ms = min(G-1, max.steps)
    print(paste0("--- Forward stepwise on response ", respind,
                 " with ", n, " obs, using ", ms, " steps"))
    results = forward_group(glX, Y, all.groups, weights, Sigma, max.steps = ms)
    Pvals[respind, 1:length(results$p.vals)] = results$p.vals
    Pvals.chi[respind, 1:length(results$chi.pvals)] = results$chi.pvals
    A.set.full = results$active.set
    A.ind = stop_last(results$p.vals)
    A.groups = A.set.full[1:A.ind]
    A.set = c()
    for (group in A.groups) {
        gname = colnames(X[,which(groups==group)])[1]
        gname = substr(gname, 1, nchar(gname)-1)
        A.set = c(A.set, gname)
    }
    A.sets[[colnames(Z)[respind]]] = A.set
}

print(A.sets)
print(sapply(A.sets, length))

filename = paste0("figs/HIV_", db, "_glint.pdf")
pdf(filename)
par(mfrow=c(2,3))
for (i in 1:nresps) {
    nonsignif = Pvals[i,(length(A.sets[[i]])+1):max.steps]
    nsignif = max.steps-length(nonsignif)
    plot(c(Pvals[i,1:nsignif], rep(NA, max.steps-nsignif)), ylim=c(0,1),
         xlab="Step", ylab="P-values", main=names(A.sets)[i])
    points(nsignif, Pvals[i,nsignif], pch=19)
    points((nsignif+1):max.steps, nonsignif, col='red')
    points((nsignif+1):max.steps, Pvals.chi[i,(nsignif+1):max.steps], col='green')
    abline(h=.1, lty='dotted')
#    clip(nsignif, nsignif, .1, 1)
#    segments(x0=nsignif+.5, y0=0.1, y1=1, lty='dotted')
}
dev.off()

texfile = paste0("tables/HIV_", db, "_glint.tex")
mat = matrix(nrow=nresps, ncol=4)
colnames(mat) = c("Drug", "n", "G", "Selected variables")
for (i in 1:nresps) {
    mat[i,1] = names(A.sets)[i]
    mat[i,2] = ns[i]
    mat[i,3] = Gs[i]
    mat[i,4] = paste(A.sets[[i]], collapse = " ")
}
caption="Variables chosen using \\textit{last} stopping rule"
tex_table(texfile, mat, caption)
