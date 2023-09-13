summary.FitDiff <- function(object, fit.measures = "default", nd = 3, tag = "\u2020") {

  if (nrow(object@nested) > 0L) {
    cat("Nested Model Comparison-------------------------------------------------------------------",
        "\n","\n")
    test.statistics <- object@nested
    if (object@model.class == "lavaan") {
      print(test.statistics, nd = nd)
    } else {
      class(test.statistics) <- c("lavaan.data.frame","data.frame")
      stats::printCoefmat(test.statistics, P.values = TRUE, has.Pvalue = TRUE)
    }
    cat("\n")
  }


  noFit <- ncol(object@fit) == 1L && names(object@fit)[1] == "df"
  if (!noFit) {
    if (is.null(fit.measures)) fit.measures <- colnames(object@fit)
    if ("all" %in% fit.measures) fit.measures <- colnames(object@fit)
    if (length(fit.measures) == 1 && fit.measures == "default") {
      ## robust or scaled test statistics?
      if (is.null(object@fit$cfi.scaled)) {
        fit.measures <- c("chisq","df","pvalue","rmsea","cfi","srmr")
      } else if (all(!is.na(object@fit$cfi.robust)) && !is.null(object@fit$cfi.robust)) {
        fit.measures <- c("chisq.scaled","df.scaled","pvalue.scaled",
                          "rmsea.robust","cfi.robust","srmr")
      } else {
        fit.measures <- c("chisq.scaled","df.scaled","pvalue.scaled",
                          "rmsea.scaled","cfi.scaled","srmr")
      }

      if ("aic" %in% colnames(object@fit)) {
        #fit.measures <- c(fit.measures, "aic", "bic")
        fit.measures <- c(fit.measures)
      }
    }


    cat("Model Fit Indices ------------------------------------------------------------------------",
        "\n","\n")
    ## this is the object to return (numeric, no printed daggers)
    fit.indices <- object@fit[ , fit.measures , drop = FALSE]

    ## print with daggers marking each fit index's preferred model
    ## (turns "numeric" vectors into "character")
    badness <- grepl(pattern = c("chisq|rmsea|ic|rmr|ecvi|fmin|hqc"),
                     x = colnames(fit.indices))
    goodness <- grepl(pattern = c("cfi|tli|rfi|nfi|ifi|rni|cn|gfi|mfi|Hat"),
                      x = colnames(fit.indices))
    minvalue <- badness & !goodness
    minvalue[!badness & !goodness] <- NA
    fit.integer <- grepl(pattern = c("df|npar|ntotal"),
                         x = colnames(fit.indices))
    suppressWarnings(fitTab <- as.data.frame(mapply(tagCharacter, nd = nd,
                                                    char = tag,
                                                    vec = fit.indices,
                                                    minvalue = minvalue,
                                                    print_integer = fit.integer),
                                             stringsAsFactors = FALSE))
    rownames(fitTab) <- object@name
    colnames(fitTab) <- colnames(fit.indices)
    class(fitTab) <- c("lavaan.data.frame","data.frame")
    print(fitTab, nd = nd)
    cat("\n")


    if (nrow(object@nested) > 0L) {
      fit.diff.measures <- fit.measures[!grepl(pattern = "chisq|pvalue|ntotal",
                                               x = fit.measures)]
      cat("Differences in Fit Indices ---------------------------------------------------------------",
          "\n","\n")
      fit.diff <- object@fit.diff[ , fit.diff.measures, drop = FALSE]
      class(fit.diff) <- c("lavaan.data.frame","data.frame")
      print(fit.diff, nd = nd)
      cat("\n")
    }
  }


  invisible(object)
}
