#db = "NRTI"
db = "PI"

binary.encoding = FALSE

max.steps = 40

setwd('..')
source('fwd_step.R')
source('generate_data.R')
source('stop_rules.R')
source('tex_table.R')


if (db == "PI") {
  data = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt", sep = "\t", header = TRUE)
  resps = 4:10
} else {
  data = read.table("http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/NRTI_DATA.txt", sep='\t', header=TRUE)
  resps = 4:9
}

nresps = length(resps)


# Clean data, restrict to cleaned subset
Z = data[,resps]

clean_col = function(column) {
    outcol = as.character(column)
    outcol[column == '.'] = '-'
    return(as.factor(outcol))
}

cnames = c()
for (j in (max(resps)+1):ncol(data)) {
    newcol = clean_col(data[,j])
    if ((length(unique(newcol)) > 1)) { 
        cnames = c(cnames, names(data)[j])
        Z = cbind(Z, newcol)
    }
}
colnames(Z) = c(names(data)[resps], cnames)

### Z should not be modified after this point ###

# Perform analyses for each possible response

A.sets = list()
ns = c()
ms = c()
Gs = c()
Pvals = matrix(NA, nrow=nresps, ncol=max.steps)
Pvals.chi = Pvals
for (respind in 1:nresps) {
    # Restrict to non-missing subset
    subrows = which(!is.na(Z[,respind]))
    Y = log(Z[subrows,respind])
    Y = Y - mean(Y)
    n = length(Y)
    # Form design matrix using non-missing subset
    X = matrix(NA, nrow=n)
    X.cat = X
    groups = c()
    g = 1
    cnames = c()
    for (j in (nresps+1):ncol(Z)) {
      newcol = Z[subrows, j]
      if (length(unique(newcol)) > 1) {
        cnames = c(cnames, names(Z)[j])
        submat = model.matrix(~ newcol - 1)
        if (binary.encoding) {
          if (ncol(submat) > 2) {
            submat = cbind(submat[,1], rowSums(submat[,2:ncol(submat)]))
          } 
          colnames(submat) = paste0(names(Z)[j], c("-", "*"))
        } else {
          colnames(submat) = paste0(names(Z)[j], levels(Z[,j]))
        }
        groups = c(groups, rep(g, ncol(submat)))
        g = g + 1
        X = cbind(X, submat)
        X.cat = cbind.data.frame(X.cat, as.factor(newcol))
      } 
    }
    X = X[ , -1]
    X.cat = X.cat[,-1]
    colnames(X.cat) = cnames
    G = max(groups)
    Gs = c(Gs, G)
    ns = c(ns, n)
    weights = sqrt(rle(groups)$lengths)
    # Estimate sigma^2 from full model
    full.model = lm(Y ~ . - 1, data=data.frame(X.cat), na.action = na.omit)
    s2 = sum(full.model$residuals^2)
    s2 = s2/full.model$df.residual
    Sigma = diag(rep(s2,n))
    # Perform forward stepwise
    m.s = min(c(G-1, max.steps, floor(n/10)))
    ms = c(ms, m.s)
    print(paste0("--- Forward stepwise on response ", respind,
                 " with ", n, " obs, using ", m.s, " steps"))
    results = forward_group(X, Y, groups, weights, Sigma, max.steps = m.s)
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
    if (length(A.set) > 0) {
      A.sets[[colnames(Z)[respind]]] = A.set
    } else {
      A.sets[[colnames(Z)[respind]]] = c(0)
    }
}

print(A.sets)
print(sapply(A.sets, length))

filename = paste0("figs/HIV_", db, "_plot.pdf")
pdf(filename)
par(mfrow=c(2,3))
for (i in 1:nresps) {
    max.steps = ms[i]
    na.vals = is.na(Pvals[i,])
    if (any(na.vals)) {
        max.steps = min(which(na.vals)) - 1
    }
    alen = length(A.sets[[i]])
    if (alen == 1) {
      if (A.sets[[i]] == 0) {
        alen = 0
      }
    }
    if (alen < max.steps) {
        # Active set is not maximum size
        nonsignif = Pvals[i,(alen+1):max.steps]
        nsignif = max.steps-length(nonsignif)
        k = nsignif
        yvals = c(Pvals[i,1:nsignif], rep(NA, max.steps-nsignif))
    } else {
        nonsignif = NULL
        nsignif = max.steps
        k = nsignif
        yvals = Pvals[i, 1:nsignif]
    }
    
    xvals = 1:length(yvals)
    plot(xvals, yvals, ylim=c(0,1), xlab="Step", ylab="P-values", main=names(A.sets)[i])
    rect(k+.5, 0, max.steps+.8, 1, col = "gray94", border = NA)
    points(xvals, yvals)
    
    if (length(nonsignif) > 0) {
        points((nsignif+1):max.steps, nonsignif)
    }
    points(1:max.steps, Pvals.chi[i,1:max.steps], col='green', pch='.', cex=3)
    abline(h=.1, lty='dotted')
}
dev.off()

texfile = paste0("tables/HIV_", db, ".tex")
mat = matrix(nrow=nresps, ncol=4)
colnames(mat) = c("Drug", "n", "G", "Selected variables")
for (i in 1:nresps) {
    mat[i,1] = names(A.sets)[i]
    mat[i,2] = ns[i]
    mat[i,3] = Gs[i]
    mat[i,4] = paste(A.sets[[i]], collapse = " ")
}
caption="Variables chosen using \\textit{last} stopping rule in forward stepwise"
tex_table(texfile, mat, caption)

## filename = paste0("figs/HIV_", db, ".pdf")
## pdf(filename)
## par(mfrow=c(2,3))
## for (i in 1:nresps) {
##     nonsignif = Pvals[i,(length(A.sets[[i]])+1):max.steps]
##     nsignif = max.steps-length(nonsignif)
##     plot(c(Pvals[i,1:nsignif], rep(NA, max.steps-nsignif)), ylim=c(0,1),
##          xlab="Step", ylab="P-values", main=names(A.sets)[i])
##     points(nsignif, Pvals[i,nsignif], pch=19)
##     points((nsignif+1):max.steps, nonsignif, col='red')
##     points((nsignif+1):max.steps, Pvals.chi[i,(nsignif+1):max.steps], col='green')
##     abline(h=.1, lty='dotted')
## #    clip(nsignif, nsignif, .1, 1)
## #    segments(x0=nsignif+.5, y0=0.1, y1=1, lty='dotted')
## }
## dev.off()

## texfile = paste0("tables/HIV_", db, ".tex")
## mat = matrix(nrow=nresps, ncol=4)
## colnames(mat) = c("Drug", "n", "G", "Selected variables")
## for (i in 1:nresps) {
##     mat[i,1] = names(A.sets)[i]
##     mat[i,2] = ns[i]
##     mat[i,3] = Gs[i]
##     mat[i,4] = paste(A.sets[[i]], collapse = " ")
## }
## caption="Variables chosen using \\textit{last} stopping rule"
## tex_table(texfile, mat, caption)
