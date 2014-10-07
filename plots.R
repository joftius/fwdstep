library(ggplot2)

step_plot = function(TrueStep, null.p, signal.p, chi.p, k, n, p, g, ugsizes, max.steps, upper.coeff, lower.coeff, max.beta, min.beta, fwd.power, design, filename,  main.append = "") {

    gmaxmin = paste0(c(max(ugsizes), min(ugsizes)), collapse="/")
    #filename = paste0('figs/', design, "_n", n, '_p', p, '_g', g,
    #    '_k', k, '_lower', lower.coeff, '_upper', upper.coeff)
    #filename = paste0(filename, fn.append)
    #filename = paste0(gsub(".", "pt", filename, fixed=TRUE), ".eps")
    filename = paste0("art/", filename, ".eps")

    ## bar.quantiles <- c(.25, .75)
    ## point.quantiles <- c(.05, .95)
    ## null.Pvals <- colMeans(null.p)
    ## null.Pvals.bar <- apply(null.p, 2, function(col) quantile(col, probs = bar.quantiles))
    ## null.Pvals.point <- apply(null.p, 2, function(col) quantile(col, probs = point.quantiles))

    ## Pvals = colMeans(signal.p)
    ## Pvals.bar <- apply(signal.p, 2, function(col) quantile(col, probs = bar.quantiles))
    ## Pvals.point <- apply(signal.p, 2, function(col) quantile(col, probs = point.quantiles))

    ## Chivals = colMeans(chi.p)
    ## Chivals.bar <- apply(chi.p, 2, function(col) quantile(col, probs = bar.quantiles))
    ## Chivals.point <- apply(chi.p, 2, function(col) quantile(col, probs = point.quantiles))

#    rownames(signal.p) = paste("signal", 1:nrow(signal.p))
    xax = 1:max.steps
    nxax = xax
    cxax = xax + .1
    plot.main = paste0("n:", n, ", p:", p, ", g:", g)
    if (g != p) {
      plot.main = paste0(plot.main, "(", gmaxmin, ")")
    }
    plot.main = paste0(plot.main, " k:", k,
                        ", beta:", upper.coeff, "/", lower.coeff,
                        "(", round(max.beta,1), "/", round(min.beta, 1),
                        "), k-TPR:", round(fwd.power, 2), main.append)

    setEPS()
    postscript(filename)

    plot(xax, TrueStep, type = "b", main = plot.main, xlab = "Step",
         ylab = "", ylim = c(0,1),
         xlim = c(min(xax) - .2, max(nxax) + .2),
         lwd=2, lty = "dashed")

    rect(k+.5, -.09, max.steps+.8, 1.09, col = "gray97", border = NA)
    lines(xax, TrueStep, type = "b", lwd=2, lty = "dashed")

#    boxplot(null.p, at=nxax, col="orangered", pch=20, add=T, boxwex=.1)
#    arrows(nxax, null.Pvals.bar[1, ], nxax, null.Pvals.bar[2, ],
#           code = 3, angle = 90, length = 0, col = "orangered")
#    points(nxax, null.Pvals.point[1, ], col = "orangered", pch = 24, cex = .5)
#    points(nxax, null.Pvals.point[2, ], col = "orangered", pch = 25, cex = .5)

    boxplot(signal.p, at=(xax - .1), col="dimgray", add=T, boxwex=.1, xaxt="n", outcex = .5)
#    arrows(xax, Pvals.bar[1, ], xax, Pvals.bar[2, ],
#           code = 3, angle = 90, length = 0, col = "blue")
#    points(xax, Pvals.point[1, ], col = "blue", pch = 24, cex = .5)
#    points(xax, Pvals.point[2, ], col = "blue", pch = 25, cex = .5)

    boxplot(chi.p, at=cxax, col="lightgray", add=T, boxwex=.1, xaxt="n", outpch=".", outcex = .5)
#    arrows(cxax, Chivals.bar[1, ], cxax, Chivals.bar[2, ],
#           code = 3, angle = 90, length = 0, col = "green")
#    points(cxax, Chivals.point[1, ], col = "green", pch = 24, cex = .5)
#    points(cxax, Chivals.point[2, ], col = "green", pch = 25, cex = .5)
    abline(h = .1)
    #axis(1, at=k, tcl=.5)
    #axis(3, at=k, tcl=.5, labels=FALSE)
    abline(v = xax[1:(length(xax)-1)] + .5)
    dev.off()

}

# Input matrix of Tchi and chisq p-values from global null simulations
qq_plot <- function(pvals, main = "", ...) {
  pvals_tchi <- sort(pvals[, 1])
  pvals_MC <- sort(pvals[, 2])
  pvals_chisq <- sort(pvals[, 3])
  n <- nrow(pvals)
  if (ncol(pvals) > 3) {
    FDR <- mean(pvals[, 4])
    main <- paste0(main, ", FDR: ", round(FDR, 2))
  }
  xrange <- seq(from = 0, to = 1, length.out = n)  
  df <- data.frame(pval = c(pvals_tchi, pvals_MC, pvals_chisq),
                   type = c(rep("tchi", n), rep("MC", n),
                       rep("chisq", n)), x = rep(xrange, 3))
  p <-
  ggplot(df, aes(x = x, y = pval)) + ggtitle(main) +
           geom_line(aes(linetype = type), size = 1) +
           geom_abline(intercept = 0, slope = 1) +
               scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("") + ylab("") + theme_bw() +
  theme(legend.position = "none", text=element_text(family = "Times"))
  return(p)
}

qq_plots <- function(df) {
  n <- nrow(df)/12
  xrange <- seq(from = 0, to = 1, length.out = n)
  if (!is.null(df$fdr)) {
    df$text <- paste0(df$text, ", FDR: ", df$fdr)
  }
  labels <- unique(cbind(df$order, df$text))
  labeltext <- labels[,2]
  labelorder <- as.numeric(labels[,1])  
  df$text <- factor(df$text, levels = labeltext[labelorder])
  df <- df[order(df$type, df$text, df$pvals),]
  df$x <- rep(xrange, 12)
  
  pl <- ggplot(df, aes(x = x, y = pvals)) +
      geom_line(aes(linetype = type), size = 1) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap( ~ text, ncol = 2) +
      scale_linetype_manual(
        values = c("dotted", "dashed", "solid")) +
      xlim(c(0,1)) + ylim(c(0,1)) +
      xlab("Uniform quantiles") +
      ylab("Observed quantiles") +
      theme_bw() +
      theme(legend.position = "none",
            text=element_text(family = "Times"),
            strip.text.x = element_text(size = 14),
            axis.text = element_text(size = 12),
            strip.background = element_rect(fill = "white"))
  return(pl)
}
