

selection_criteria <- function(model, truth) {
  if ((length(model) == 1) && (model == 0)) {
    #out <- c(0, 0, 0)
    return(list(R = 0, TPR = 0, PPV = 0))
  } else {
    TP <- sum(model %in% truth)
    P <- length(truth)
    R <- length(model)
    return(list(R = R, TPR = TP/P, PPV = TP/R))
    #out <- c(R, TP/P, TP/R)
  }
#  names(out) <- c("R", "TPR", "PPV")
#  return(out)
}

