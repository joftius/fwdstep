db = "PI"
#db = "NRTI"

binary.encoding = TRUE

max.steps = 40

setwd('..')
source('fwd_step.R')
source('generate_data.R')
source('glint/generate_glint.R')
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
    if (length(unique(newcol)) > 1) {
        cnames = c(cnames, names(data)[j])
        Z = cbind(Z, newcol)
    }
}
colnames(Z) = c(names(data)[resps], cnames)
# Z should not be modified after this point

# Perform analyses for each possible response

A.sets = list()
ns = c()
Gs = c()
ps = c()
ms = c()
omits = c()
Pvals = matrix(NA, nrow=nresps, ncol=max.steps)
Pvals.chi = Pvals
for (respind in 1:nresps) {
    subrows = which(!is.na(Z[,respind]))
    Y = log(Z[subrows,respind])
    Y = Y - mean(Y)
    n = length(Y)
    X = matrix(NA, nrow=n)
    X.cat = X
    groups = c()
    g = 1
    cnames = c()
    for (j in (nresps+1):ncol(Z)) {
        newcol = Z[subrows, j]
        if (length(unique(newcol)) > 1) {
            cnames = c(cnames, names(Z)[j])
            submat = model.matrix(~ Z[subrows,j] - 1)
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
    X.cat = X.cat[, -1]
    colnames(X.cat) = cnames
    full.model = lm(Y ~ . - 1, data=data.frame(X.cat), na.action = na.omit)
    s2 = sum(full.model$residuals^2)
    s2 = s2/full.model$df.residual
    Sigma = diag(rep(s2,n))
    
    gld = generate_glint(X, groups, cat.groups = unique(groups))
    glX = gld$X
    all.groups = gld$all.groups
    special.groups = gld$special.groups
    default.groups = gld$default.groups
    omits = c(omits, gld$omitted)
    
    #glX = scale(glX, center=TRUE, scale=FALSE)
    glX = frob_normalize(glX, all.groups)
    G = max(all.groups)
    Gs = c(Gs, G)
    ns = c(ns, n)
    ps = c(ps, length(all.groups))
    #weights = sqrt(rle(groups)$lengths)
    weights = rep(1, max(all.groups))

    # Perform forward stepwise
    m.steps = min(c(G-1, max.steps, floor(n/10)))
    ms = c(ms, m.steps)
    print(paste0("--- Response ", respind, "(", round(s2, 2), ")",
                 " with ", n, " obs and ", G, " groups, using ",
                 m.steps, " steps"))
    results = forward_group(glX, Y, all.groups, weights, Sigma, max.steps = m.steps)
    Pvals[respind, 1:length(results$p.vals)] = results$p.vals
    Pvals.chi[respind, 1:length(results$chi.pvals)] = results$chi.pvals
    A.set.full = results$active.set
    A.ind = stop_last(results$p.vals)
    A.groups = A.set.full[1:A.ind]
    A.set = c()
    for (group in A.groups) {
      if (group <= max(groups)) {
        gname = colnames(X[,which(groups==group)])[1]
        gname = substr(gname, 1, nchar(gname)-1)
      } else {
        mgs = glint_default_group_of(group, default.groups)
        gname1 = colnames(X[,which(groups==mgs[1])])[1]
        gname1 = substr(gname1, 1, nchar(gname1)-1)
        gname2 = colnames(X[,which(groups==mgs[2])])[1]
        gname2 = substr(gname2, 1, nchar(gname2)-1)
        gname = paste0(gname1, "*", gname2)
      }
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

filename = paste0("figs/HIV_", db, "_glint.pdf")
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

texfile = paste0("tables/HIV_", db, "_glint.tex")
mat = matrix(nrow=nresps, ncol=4)
colnames(mat) = c("Drug", "n", "G", "Selected variables")
for (i in 1:nresps) {
    mat[i,1] = names(A.sets)[i]
    mat[i,2] = ns[i]
    mat[i,3] = Gs[i]
    mat[i,4] = paste(A.sets[[i]], collapse = " ")
}
caption="Variables chosen using \\textit{last} stopping rule in \\textsc{Glinternet}"
tex_table(texfile, mat, caption)
