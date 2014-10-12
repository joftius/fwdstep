library(MASS)
library(Matrix)
library(devtools)
library(parallel)
load_all("package/fstep/")
source("plots.R") # for qq_plot

# Load data
library(lars)
data(diabetes)
x <- diabetes$x
labels <- 1:ncol(x)
index <- labels
X <- scale_groups(x, index, scale = FALSE)
Xnorm <- scale_groups(X, index)
Y <- diabetes$y
Y <- Y - mean(Y)

# Scale Y
bigmodel <- lm(Y ~ Xnorm)
Y <-  Y/sqrt(mean(resid(bigmodel)^2))

# Run forward stepwise
fit <- fstep(x=Xnorm, y=Y, Sigma=1, steps = 9)

# Create table of output
pvals <- data.frame(round(rbind(fit$p.value, fit$log$chisq), 3))
colnames(pvals) <- colnames(X)[fit$variable]
rownames(pvals) <- c("Tchi", "chisq")

xtable(pvals, caption = "Forward stepwise results on diabetes dataset. Chi-squared $p$-values are computed using the drop in residual sum of squares. Variables are listed in the order of inclusion by forward stepwise.", digits = 3)

