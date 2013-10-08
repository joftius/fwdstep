
source('pred_and_estim.R')

main = function(num.nonzero = 10) {
  n = 100
  for (n in c(100, 500)) {
    niter = 200
    stop.rule = "first"
    X = diag(rep(1, n))
    X.test = X
    groups = 1:n
    beta = c(rep(1, num.nonzero), rep(0, n - num.nonzero))
  
    for (signal in c(1, 3.5)) {
      outmat = matrix(ncol = 3)
      for (iter in 1:niter) {
        train.noise = rnorm(n)
        test.noise = rnorm(n)
        Y = beta * signal + train.noise
        Y.test = beta * signal + test.noise
        orderY = sort(abs(Y), decreasing = TRUE)
        p.list = 2*pnorm(orderY, lower.tail = FALSE)
        active.set = rev(order(abs(Y)))
        stats = pred_est_stats(p.list, active.set, X, Y, groups, signal * beta, X.test, Y.test, stop.rule)
        outmat = rbind(outmat, unlist(stats))
      }
      outmat = outmat[1 + 1:niter, ]
      print(c(colMeans(outmat), apply(outmat, 2, sd)))
    }
  }
}

main()

main(50)
