library(MASS)
library(Matrix)
library(devtools)
library(parallel)
load_all("package/fstep/")
source("plots.R") # for qq_plot
source("generate_data.R") # for beta_staircase


set.seed(1)

nsim <- 500
high <- 15
low <- 12

# ----------------------------------------------------------

library(flare)
data(eyedata)
labels <- 1:ncol(x)
index <- labels
X <- scale_groups(x, index, scale = FALSE)
Xnorm <- scale_groups(X, index)

# ----------------------------------------------------------

n <- nrow(X)
p <- length(index)
G <- length(labels)

instance <- function(n, p, index, labels, k, steps) {

  if (k > 0) {
    mult          <- sqrt(2*log(p)/n)
    upper         <- high * mult
    lower         <- low * mult
    bd            <- beta_staircase(groups = index, num.nonzero = k, upper = upper, lower = lower)
    beta          <- bd$beta
    truth         <- bd$true.active
    Y     <- X %*% beta + rnorm(n)

  } else {
    Y <- rnorm(n)
  }

  Y     <- Y - mean(Y)
  fit   <- fstep(x = Xnorm, y = Y, index, Sigma = 1, steps = steps)
  active <- fit$variable
  inactive <- setdiff(labels, active)
  inactive.inds <- which(!(index %in% active))


  if (k > 0) {
    FP <- length(setdiff(active, truth))
    FDR <- FP/steps
  }

  chisq.pval <- fit$log$chisq[steps]

  Xnorm.adj <- Xnorm
  cand.labels <- c(inactive, active[steps])
  L <- fit$log$L[steps]
  if (steps > 1) {
    adj.inds <- which(index %in% active[1:(steps-1)])
    cand.inds <- which(index %in% cand.labels)
    Xnorm.adj[, cand.inds] <- resid(lm(Xnorm[, cand.inds] ~ Xnorm[, adj.inds] -1))
    lastY <- resid(lm(Y ~ Xnorm[, adj.inds] - 1))
    L <- t(Xnorm[, which(index == active[steps])]) %*% lastY
    L <- sqrt(sum(L^2))
  }

  mc.pval <- MC_pvalue(L, Xnorm.adj, index, cand.labels)

  output <- data.frame(
      pvals = c(fit$p.value[steps], mc.pval, chisq.pval),
      type = c("tchi", "mc", "chisq"))
  if (k > 1) {
    output$fdr <- FDR
  } else {
    output$fdr <- 1
  }  
  return(output)
}

simulation <- function(nsim, n, p, index, labels, k, steps){
  results <- replicate(nsim, instance(n, p, index, labels, k, steps), simplify = FALSE)
  output <- do.call(rbind, results)
  output$fdr <- round(mean(output$fdr), 3)  
  return(output)
}


four_four <- simulation(nsim, n, p, index, labels, k = 0, steps = 1)
four_four$text <- "Null, step 1"
four_four$order <- 1

four_five <- simulation(nsim, n, p, index, labels, k = 4, steps = 4)
four_five$text <- "4-sparse, step 4"
four_five$order <- 2

four_six <- simulation(nsim, n, p, index, labels, k = 4, steps = 5)
four_six$text <- "4-sparse, step 5"
four_six$order <- 3

four_seven <- simulation(nsim, n, p, index, labels, k = 4, steps = 6)
four_seven$text <- "4-sparse, step 6"
four_seven$order <- 4


df <- data.frame(rbind(four_four, four_five, four_six, four_seven))

setEPS()
postscript(paste0("art/nonnull_G", G, "_high", high, "_data.eps"), fonts = "Times")

qq_plots(df)

dev.off()


