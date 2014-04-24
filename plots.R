

step_plot = function(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper.coeff, lower.coeff, max.beta, min.beta, fwd.power, design, filename,  main.append = "") {

    gmaxmin = paste0(c(max(ugsizes), min(ugsizes)), collapse="/")
    #filename = paste0('figs/', design, "_n", n, '_p', p, '_g', g,
    #    '_k', k, '_lower', lower.coeff, '_upper', upper.coeff)
    #filename = paste0(filename, fn.append)
    #filename = paste0(gsub(".", "pt", filename, fixed=TRUE), ".pdf")
    filename = paste0("figs/", filename, ".pdf")

    bar.quantiles <- c(.25, .75)
    point.quantiles <- c(.05, .95)
    null.Pvals <- colMeans(null.p)
    null.Pvals.bar <- apply(null.p, 2, function(col) quantile(col, probs = bar.quantiles))
    null.Pvals.point <- apply(null.p, 2, function(col) quantile(col, probs = point.quantiles))

    Pvals = colMeans(signal.p)
    Pvals.bar <- apply(signal.p, 2, function(col) quantile(col, probs = bar.quantiles))
    Pvals.point <- apply(signal.p, 2, function(col) quantile(col, probs = point.quantiles))

    Chivals = colMeans(chi.p)
    Chivals.bar <- apply(chi.p, 2, function(col) quantile(col, probs = bar.quantiles))
    Chivals.point <- apply(chi.p, 2, function(col) quantile(col, probs = point.quantiles))

    nxax <- 1:max.steps - 0.1
    xax <- nxax - 0.1
    cxax = nxax - 0.2
    plot.main <- paste0("n:", n, ", p:", p, ", g:", g)
    if (g != p) {
      plot.main = paste0(plot.main, "(", gmaxmin, ")")
    }
    plot.main = paste0(plot.main, " k:", k,
                        ", beta:", upper.coeff, "/", lower.coeff,
                        "(", round(max.beta,1), "/", round(min.beta, 1),
                        "), k-SP:", round(fwd.power, 2), main.append)

    pdf(filename)

    plot(xax, TrueStep, type = "b", main = plot.main, xlab = "Step",
         ylab = "", ylim = c(-.01,1.01),
         xlim = c(min(xax) - .2, max(nxax) + .2),
         lwd=2, lty = "dashed")
    rect(k+.5, -.05, max.steps+.8, 1.05, col = "gray94", border = NA)
    lines(xax, TrueStep, type = "b", lwd=2, lty = "dashed")

    points(nxax, null.Pvals, col="orangered", pch=20)
    arrows(nxax, null.Pvals.bar[1, ], nxax, null.Pvals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "orangered")
    points(nxax, null.Pvals.point[1, ], col = "orangered", pch = 24, cex = .5)
    points(nxax, null.Pvals.point[2, ], col = "orangered", pch = 25, cex = .5)

    points(xax, Pvals, col="blue", pch=15)
    arrows(xax, Pvals.bar[1, ], xax, Pvals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "blue")
    points(xax, Pvals.point[1, ], col = "blue", pch = 24, cex = .5)
    points(xax, Pvals.point[2, ], col = "blue", pch = 25, cex = .5)

    points(cxax, Chivals, col="green", pch=18)
    arrows(cxax, Chivals.bar[1, ], cxax, Chivals.bar[2, ],
           code = 3, angle = 90, length = 0, col = "green")
    points(cxax, Chivals.point[1, ], col = "green", pch = 24, cex = .5)
    points(cxax, Chivals.point[2, ], col = "green", pch = 25, cex = .5)
    abline(h = .1)
    #axis(1, at=k, tcl=.5)
    #axis(3, at=k, tcl=.5, labels=FALSE)

    dev.off()
}
