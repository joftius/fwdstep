library(MASS)
library(Matrix)
library(devtools)
library(ggplot2)
library(gridExtra)
load_all("package/fstep/")
source("plots.R")
source("selection.R")
source("generate_data.R")

set.seed(1)

nsim <- 500
n <- 100
p <- 200
G <- p
maxk <- 8
index <- 1:p
labels <- unique(index)
high <- 2
low <- 1.5

instance <- function(n, p, maxk, index, labels) {
  X     <- matrix(rnorm(n*p), nrow = n)
  Xnorm <- scale_groups(X, index)
  G     <- length(unique(labels))
  noise <- rnorm(n)
  klist <- 2*(1:maxk)
  out   <- lapply(klist, function(k) k_instance(k, noise, X, Xnorm, index, labels))
  out   <- data.frame(do.call(rbind, out))
  return(out)
}

k_instance <- function(k, noise, X, Xnorm, index, labels) {
  n             <- nrow(X)
  p             <- ncol(X)
  rule_names    <- c("last", "forward", "RIC", "BIC")
  mult          <- sqrt(2*log(p)/n)
  upper         <- high * mult
  lower         <- low * mult
  bd            <- beta_staircase(index, k, upper, lower)
  beta          <- bd$beta
  truth         <- bd$true.active
  Y             <- X %*% beta + noise
  fit           <- fstep(x=Xnorm, y=Y, Sigma=1, steps = floor(n/2))
  models        <- model_select(fit)
  out           <- lapply(models, function(model) selection_criteria(model, truth))
  out           <- data.frame(do.call(rbind, out))
  out$rule      <- factor(rule_names)
  out$k         <- k
  rownames(out) <- NULL
  return(out)
}

simulation <- function(nsim, n, p, maxk, index, labels){
  results <- replicate(nsim, instance(n, p, maxk, index, labels), simplify = FALSE)
  df <- data.frame(do.call(rbind, results))
  means <- aggregate(df[,1:3], by = list(df$rule, df$k),
                     function(coercing) mean(as.numeric(coercing)))
  sdevs <- aggregate(df[,1:3], by = list(df$rule, df$k),
                        function(coercing) sd(as.numeric(coercing)))
  colnames(means)[1:2] <- c("rule", "k")
  colnames(sdevs)[1:2] <- c("rule", "k")
  return(list(means = means, deviations = sdevs))
}

data <- simulation(nsim, n, p, maxk, index, labels)

dm <- data$means
ds <- data$deviations

p1 <- ggplot(dm, aes(x=k, y=TPR, linetype=rule)) + ylim(0,1) +
    xlab("Sparsity") + geom_line() + theme_bw() +
    scale_linetype_manual(values = c("dotted", "dashed", "dotdash", "solid")) + theme(legend.position = "none")
p2 <- ggplot(dm, aes(x=k, y=PPV, linetype=rule)) + ylim(0,1) +
    xlab("Sparsity") + geom_line() + theme_bw() +
    scale_linetype_manual(values = c("dotted", "dashed", "dotdash", "solid")) + ylab("1 - FDR") +
    theme(legend.position = "none",
          text=element_text(family = "Times"))

setEPS()
postscript(paste0("art/selection_G", G, "_high", high, ".eps"), width = 20, height = 10, paper = "special", fonts = "Times")

grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size="last"))

dev.off()

write.csv(dm, "data/selection_means_high.csv")
write.csv(ds, "data/selection_sdevs_high.csv")

