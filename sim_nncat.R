library(MASS)
library(Matrix)
library(devtools)
library(parallel)
load_all("package/fstep/")
source("plots.R") # for qq_plot
source("generate_data.R") # for beta_staircase

set.seed(1)

nsim <- 500
n <- 100
G <- 100
high <- 4
low <- 2
labels <- 1:G
index <- sort(rep(labels, 2))
p <- length(index)


instance <- function(n, p, index, labels, k, steps) {

  X <- matrix(nrow=n, ncol=p)
  betaprior <- rbeta(G, sqrt(G), sqrt(G))
  for (i in 1:G) {
      X[, 2*i-1] <- sample(c(0,1), n, replace=TRUE, prob=c(betaprior[i], 1 - betaprior[i]))
      X[, 2*i] <- 1 - X[, 2*i-1]
  }
  Xnorm <- X/sqrt(n)

  if (k > 0) {
    mult          <- sqrt(2*log(p)/n)
    upper         <- high * mult
    lower         <- low * mult
    bd            <- beta_staircase(groups = index, num.nonzero = k, upper = upper, lower = lower, cat.groups = labels)
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

five_four <- simulation(nsim, n, p, index, labels, k = 5, steps = 4)
five_four$text <- "5-sparse, step 4"
five_four$order <- 1

five_five <- simulation(nsim, n, p, index, labels, k = 5, steps = 5)
five_five$text <- "5-sparse, step 5"
five_five$order <- 2

five_six <- simulation(nsim, n, p, index, labels, k = 5, steps = 6)
five_six$text <- "5-sparse, step 6"
five_six$order <- 3

five_seven <- simulation(nsim, n, p, index, labels, k = 5, steps = 7)
five_seven$text <- "5-sparse, step 7"
five_seven$order <- 4

df <- data.frame(rbind(five_four, five_five, five_six, five_seven))

setEPS()
postscript(paste0("art/nonnull_G", G, "_high", high, "_binary.eps"), fonts = "Times")

qq_plots(df)

dev.off()

