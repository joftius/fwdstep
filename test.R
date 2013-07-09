###############################################
# Run a battery of tests to check development #
###############################################

source("generate_data.R")
source("forward_step.R")

n <- 100
g <- 20
k <- 5

data <- simulate_fixed(n, g, k)

forward_group(data$X, data$Y, data$groups, weights="default")

