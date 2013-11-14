
source('fstep_minimal.R')


n = 10
groups = c(1,1,2,2,3,4,5,6,7,8,9,10,10,10,11,11,12,12,13,13,
  14,15,16,17,17,18,18,19,19,20,20,21,22,23,24,25,26,26,26,26,26)
p = length(groups)
beta = as.numeric(groups == 2) + as.numeric(groups == 6)
beta = sqrt(2*log(p))*beta
max.steps = 2

for (i in 1:10) {
  X = matrix(rnorm(n*p), nrow=n)
  Y = X %*% beta
  noise = rnorm(n)
  print(fstep_fit(X, Y + noise, groups, max.steps))
}

