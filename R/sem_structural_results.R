sem_structural_results <- function(model, nd = 3) {
  indic <- which(inspect(model, what = "std")$beta != 0,
                 arr.ind = TRUE, useNames = TRUE)

  result <- as.data.frame(cbind(
    colnames(inspect(model, what = "std")$beta)[indic[, 2]],
    colnames(inspect(model, what = "std")$beta)[indic[, 1]],
    round(inspect(model, what = "std")$beta[indic], digits = 3),
    round(inspect(model, what = "se")$beta[indic], digits = 3),
    round(pnorm(abs(inspect(model, what = "coef")$beta /
                      inspect(model, what = "se")$beta), lower.tail = FALSE)[indic], digits = 3)
  )
  )
  colnames(result) <- c("outcome", "predictor", "std estimate", "se", "p-value")

  return(result)
}

