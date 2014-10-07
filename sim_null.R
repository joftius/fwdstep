library(MASS)
library(Matrix)
library(devtools)
library(ggplot2)
library(gridExtra)
load_all("package/fstep/")
source("plots.R")

set.seed(1)

nsim <- 500
n <- 100
p <- 200
G <- p
index <- 1:p
labels <- unique(index)

instance <- function(type, n, p, index, labels, wrong = 1) {
  X <- scale_groups(matrix(rnorm(n*p), nrow = n), index)
  if (type == "sigma") {
    Y <- wrong * rnorm(n)
  } else if (type == "heavy") {
    Y <- rt(n, df = 18/5)
  } else if (type == "skew") {
    f <- function(s) exp(s^2) * (exp(s^2) - 1) - 1/2
    s <- uniroot(f, c(.5, .6))$root
    Y <- rnorm(n) + rlnorm(n, sdlog = s)
  }

  Y <- Y - mean(Y)
  fit <- add1.fstep(x = X, y = Y, index = index, Sigma = 1)
  inactive <- setdiff(labels, fit$imax)

  RSSdrop <- sum(Y^2) - sum(resid(lm(Y ~ X[, which(index == fit$imax)]))^2)

  chisq.pval <- pchisq(RSSdrop, lower.tail=FALSE, df = 1)
  mc.pval <- MC_pvalue(fit$L, X, index, inactive)

  return(data.frame(
      pvals = c(fit$p.value, mc.pval, chisq.pval),
      type = c("tchi", "mc", "chisq")))
}

simulation <- function(type, nsim, n, p, index, labels, wrong){
  results <- replicate(nsim, instance(type, n, p, index, labels, wrong), simplify = FALSE)
  return(do.call(rbind, results))
}

largesigma <- simulation(type = "sigma", nsim, n, p, index, labels, wrong = 0.8)
largesigma$text <- "Sigma overstimated"
largesigma$order <- 1

smallsigma <- simulation(type = "sigma", nsim, n, p, index, labels, wrong = 1.2)
smallsigma$text <- "Sigma underestimated"
smallsigma$order <- 2

heavy <- simulation(type = "heavy", nsim, n, p, index, labels)
heavy$text <- "Heavy-tailed"
heavy$order <- 3

skewed <- simulation(type = "skew", nsim, n, p, index, labels)
skewed$text <- "Skewed"
skewed$order <-  4

df <- data.frame(rbind(largesigma, smallsigma, heavy, skewed))

setEPS()
postscript(paste0("art/global_null_G", G, ".eps"), fonts = "Times")

qq_plots(df)

dev.off()


