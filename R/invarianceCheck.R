invarianceCheck <- function(model, data, group, estimator = "MLR",
                            intercept = FALSE, missing = FALSE, display = FALSE) {

  if (estimator == "MLR") {
    if (intercept == FALSE) {
      if (missing == FALSE) {
        cfa <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                   group = group, estimator = "MLR")

        loading_invariance <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                                  group = group,
                                  group.equal = c("loadings"),
                                  estimator = "MLR")
      } else{
        cfa <- cfa(model, data, std.lv = TRUE, missing  = "fiml",
                   group = group, estimator = "MLR")

        loading_invariance <- cfa(model, data, std.lv = TRUE,
                                  missing = "fiml", group = group,
                                  group.equal = c("loadings"), estimator = "MLR")
      }

      load_check <- semTools::compareFit(Loadings_Invariant = loading_invariance,
                                         Loadings_Free = cfa)
      summary.FitDiff(load_check)

      pvalue_chisq_diff <- as.numeric(
        anova(loading_invariance, cfa)$"Pr(>Chisq)"[2])

      diff_cfi_scaled <- fitMeasures(loading_invariance)[["cfi.scaled"]] -
        fitMeasures(cfa)[["cfi.scaled"]]
      diff_rmsea_scaled <- fitMeasures(loading_invariance)[["rmsea.scaled"]] -
        fitMeasures(cfa)[["rmsea.scaled"]]
      diff_srmr_scaled <- fitMeasures(loading_invariance)[["srmr"]] -
        fitMeasures(cfa)[["srmr"]]

      # Chi-Square --------------------------------------------------------------
      cat("Loading Invariance Interpretation --------------------------------------------------------", "\n", "\n")

      if (pvalue_chisq_diff >= 0.05) {
        cat("The hypothesis of perfect loading invariance *is not* rejected according to the
      Chi-Square difference test statistics because the p-value is larger or equal to 0.05.", "\n", "\n")
      }else {
        cat("The hypothesis of perfect loading invariance *is* rejected according
           to the Chi-Square difference test statistics because the p-value is
           smaller than 0.05", "\n", "\n")
      }
      # CFI  --------------------------------------------------------------------
      if (diff_cfi_scaled <= .01) {
        cat("The hypothesis of approximate loading invariance *is not* rejected
        according to the CFI because the difference in CFI value is smaller than
        or equal to 0.01.", "\n", "\n")
      }else {
        cat("The hypothesis of approximate approximate loading invariance *is*
        rejected according to the CFI because the difference in CFI value is
        larger than 0.01.", "\n", "\n")
      }


      # RMSEA -------------------------------------------------------------------
      if (diff_rmsea_scaled <= 0.015) {
        cat("The hypothesis of approximate loading invariance *is not* rejected
       according to the RMSEA because the difference in RMSEA value is smaller
       than or equal to 0.015.", "\n", "\n")
      } else {
        cat("The hypothesis of approximate approximate loading invariance *is*
       rejected according to the RMSEA because the difference in RMSEA value is
       larger than 0.015.", "\n", "\n")
      }

      # SRMR --------------------------------------------------------------------
      if (diff_srmr_scaled < .030) {
        cat("The hypothesis of approximate loading invariance *is not* rejected
       according to the SRMR because the difference in SRMR value is smaller than
       or equal to 0.030.", "\n", "\n")
      } else {
        cat("The hypothesis of approximate approximate loading invariance *is*
       rejected according to the SRMR because the difference in SRMR value is
       larger than 0.030.", "\n", "\n")
      }

      t1 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                       as.numeric(levels(pull(data, group))[1]),"`$lambda)")
      )))

      loadings <- array(NA, dim = c(t1[1], t1[2],length(levels(pull(data, group)))))
      group_upd <- as.numeric(pull(data, group))

      for (i in min(group_upd):max(group_upd)){
        eval(parse(text = c(
          paste0("loadings[,,", i, "] <-", "inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[i]), "`$lambda")
        )))
      }

      rowlab <- eval(parse(text = c(paste0("dimnames(inspect(cfa, \"est\")$`",
                                           as.numeric(levels(pull(data, group))[1]),"`$lambda)[[1]]")
      )))

      collab <- eval(parse(text = c(paste0("dimnames(inspect(cfa, \"est\")$`",
                                           as.numeric(levels(pull(data, group))[1]),"`$lambda)[[2]]")
      )))

      dimnames(loadings) <- list(rowlab, collab, levels(pull(data, group)))

      combis <- as.matrix(combn(unique(group_upd), 2))

      loading_diff <- array(NA, dim = c(t1[1], t1[2], ncol(combis)))
      label_diff <- rep(NA, ncol(combis))

      for (j in 1:ncol(combis)){
        loading_diff[,,j] <- loadings[,,combis[1,j]] - loadings[,,combis[2,j]]
        label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                levels(pull(data, group))[combis[2,j]])
      }

      dimnames(loading_diff) <- list(rowlab, collab, label_diff)
      load_diff_abs <- abs(loading_diff)
      maxind <- which(load_diff_abs == max((load_diff_abs)), arr.ind = TRUE)

      if (display == TRUE) {
        cat("Factor loadings for each group ---------------------------------------------------------", "\n", "\n")
        print(loadings)
        cat("Difference in loadings between groups --------------------------------------------------", "\n", "\n")
        print(loading_diff)

        print(paste("The largest difference in loading is that of variable",
                    dimnames(loading_diff)[[1]][maxind[1]],"on factor",
                    dimnames(loading_diff)[[2]][maxind[2]], "for group comparison",
                    dimnames(loading_diff)[[3]][maxind[3]]))
      }

    } else {
      if (missing == FALSE) {
        cfa <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                   group = group, estimator = "MLR")

        loading_invariance <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                                  group = group,
                                  group.equal = c("loadings"), estimator = "MLR")

        intercept_invariance <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                                    group = group,
                                    group.equal = c("loadings", "intercepts"),
                                    estimator = "MLR")

      } else {
        cfa <- cfa(model, data, std.lv = TRUE, missing  = "fiml",
                   group = group, estimator = "MLR")

        loading_invariance <- cfa(model, data, std.lv = TRUE,
                                  missing = "fiml", group = group,
                                  group.equal = c("loadings"), estimator = "MLR")

        intercept_invariance <- cfa(model, data, std.lv = TRUE,
                                    missing = "fiml", group = group,
                                    group.equal = c("loadings", "intercepts"),
                                    estimator = "MLR")
      }


      pvalue_chisq_diff1 <- as.numeric(
        anova(intercept_invariance, loading_invariance, cfa)$"Pr(>Chisq)"[2])

      diff_cfi_scaled1 <- fitMeasures(loading_invariance)[["cfi.scaled"]] -
        fitMeasures(cfa)[["cfi.scaled"]]
      diff_rmsea_scaled1 <- fitMeasures(loading_invariance)[["rmsea.scaled"]] -
        fitMeasures(cfa)[["rmsea.scaled"]]
      diff_srmr_scaled1 <- fitMeasures(loading_invariance)[["srmr"]] -
        fitMeasures(cfa)[["srmr"]]

      SC1 <- sum(as.numeric(pvalue_chisq_diff1 >.05),
                 as.numeric(diff_cfi_scaled1 <= .01),
                 as.numeric(diff_rmsea_scaled1 <= 0.015),
                 as.numeric(diff_srmr_scaled1 <= 0.030))

      if (SC1 >= 3){
        int_check <- semTools::compareFit(Intercept_Invariant =
                                            intercept_invariance, Loadings_Invariant = loading_invariance,
                                          Loadings_Intercept_Free = cfa)

        summary.FitDiff(int_check)

        pvalue_chisq_diff2 <- as.numeric(
          anova(intercept_invariance, loading_invariance, cfa)$"Pr(>Chisq)"[3])

        diff_cfi_scaled2 <- fitMeasures(intercept_invariance)[["cfi.scaled"]] -
          fitMeasures(loading_invariance)[["cfi.scaled"]]
        diff_rmsea_scaled2 <- fitMeasures(intercept_invariance)[["rmsea.scaled"]] -
          fitMeasures(loading_invariance)[["rmsea.scaled"]]
        diff_srmr_scaled2 <- fitMeasures(intercept_invariance)[["srmr"]] -
          fitMeasures(loading_invariance)[["srmr"]]

        # Chi-Square ------------------------------------------------------------
        cat("Loading Invariance Interpretation --------------------------------------------------------", "\n", "\n")

        if (pvalue_chisq_diff1 >= 0.05){
          cat("The hypothesis of perfect loading invariance *is not* rejected
          according to the Chi-Square difference test statistics because the p-value
          is larger or equal to 0.05.", "\n", "\n")
        } else {
          cat("The hypothesis of perfect loading invariance *is* rejected according
          to the Chi-Square difference test statistics because the p-value is
          smaller than 0.05", "\n", "\n")
        }

        # CFI  ------------------------------------------------------------------
        if (diff_cfi_scaled1 <= .01) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
         according to the CFI because the difference in CFI value is smaller than
         or equal to 0.01.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
         rejected according to the CFI because the difference in CFI value is
         larger than 0.01.", "\n", "\n")
        }

        # RMSEA -----------------------------------------------------------------
        if(diff_rmsea_scaled1 <= 0.015) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
         according to the RMSEA because the difference in RMSEA value is smaller
         than or equal to 0.015.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
         rejected according to the RMSEA because the difference in RMSEA value is
         larger than 0.015.", "\n", "\n")
        }

        # SRMR ------------------------------------------------------------------
        if (diff_srmr_scaled1 <= .030) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
         according to the SRMR because the difference in SRMR value is smaller than
         or equal to 0.030.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
         rejected according to the SRMR because the difference in SRMR value is
         larger than 0.030.", "\n", "\n")
        }

        cat("Intercept Invariance Interpretation ------------------------------------------------------", "\n", "\n")

        if (pvalue_chisq_diff2 >= 0.05){
          cat("The hypothesis of perfect intercept invariance *is not* rejected
         according to the Chi-Square difference test statistics because the p-value
         is larger or equal to 0.05.", "\n", "\n")
        } else {
          cat("The hypothesis of perfect intercept invariance *is* rejected
         according to the Chi-Square difference test statistics because the p-value
         is smaller than 0.05", "\n", "\n")
        }

        if (diff_cfi_scaled2 <= .01) {
          cat("The hypothesis of approximate intercept invariance *is not*
         rejected according to the CFI because the difference in CFI value is
         smaller than or equal to 0.01.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate intercept invariance *is*
         rejected according to the CFI because the difference in CFI value is
         larger than 0.01.", "\n", "\n")
        }

        if (diff_rmsea_scaled2 <= 0.015) {
          cat("The hypothesis of approximate intercept invariance *is not* rejected
         according to the RMSEA because the difference in RMSEA value is smaller
         than or equal to 0.015.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate intercept invariance *is*
         rejected according to the RMSEA because the difference in RMSEA value is
         larger than 0.015.", "\n", "\n")
        }

        if (diff_srmr_scaled2 <= .030) {
          cat("The hypothesis of approximate intercept invariance *is not* rejected
         according to the SRMR because the difference in SRMR value is smaller than
         or equal to 0.030.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate intercept invariance *is*
         rejected according to the SRMR because the difference in SRMR value is
         larger than 0.030.", "\n", "\n")
        }

        t1 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                         as.numeric(levels(pull(data, group))[1]),"`$lambda)")
        )))

        loadings <- array(NA, dim = c(t1[1], t1[2],length(levels(pull(data, group)))))
        group_upd <- as.numeric(pull(data, group))

        for (i in min(group_upd):max(group_upd)){
          eval(parse(text = c(
            paste0("loadings[,,", i, "] <-", "inspect(cfa, \"est\")$`",
                   as.numeric(levels(pull(data, group))[i]), "`$lambda")
          )))
        }

        rowlab <- eval(parse(text = c(paste0("dimnames(inspect(cfa, \"est\")$`",
                                             as.numeric(levels(pull(data, group))[1]),"`$lambda)[[1]]")
        )))

        collab <- eval(parse(text = c(paste0("dimnames(inspect(cfa, \"est\")$`",
                                             as.numeric(levels(pull(data, group))[1]),"`$lambda)[[2]]")
        )))

        dimnames(loadings) <- list(rowlab, collab, levels(pull(data, group)))

        combis <- as.matrix(combn(unique(group_upd), 2))

        loading_diff <- array(NA, dim = c(t1[1], t1[2], ncol(combis)))
        label_diff <- rep(NA, ncol(combis))

        for (j in 1:ncol(combis)){
          loading_diff[,,j] <- loadings[,,combis[1,j]] - loadings[,,combis[2,j]]
          label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                  levels(pull(data, group))[combis[2,j]])
        }

        dimnames(loading_diff) <- list(rowlab, collab, label_diff)


        t2 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                         as.numeric(levels(pull(data, group))[1]),"`$nu)")
        )))

        intercepts <- array(NA, dim = (c(t2[1], length(levels(pull(data, group))))))
        group_upd <- as.numeric(pull(data, group))

        for (i in min(group_upd):max(group_upd)){
          eval(parse(text = c(
            paste0("intercepts[,", i, "] <-", "inspect(cfa, \"est\")$`",
                   as.numeric(levels(pull(data, group))[i]), "`$nu")
          )))
        }

        rowlab2 <- eval(parse(text = c(paste0("dimnames(inspect(cfa, \"est\")$`",
                                              as.numeric(levels(pull(data, group))[1]),"`$nu)[[1]]")
        )))

        dimnames(intercepts) <- list(rowlab2, levels(pull(data, group)))

        combis <- as.matrix(combn(unique(group_upd), 2))

        intercept_diff <- array(NA, dim = c(t2[1], ncol(combis)))
        label_diff <- rep(NA, ncol(combis))

        for (j in 1:ncol(combis)){
          intercept_diff[,j] <- intercepts[,combis[1,j]] - intercepts[,combis[2,j]]
          label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                  levels(pull(data, group))[combis[2,j]])
        }

        dimnames(intercept_diff) <- list(rowlab2, label_diff)

        load_diff_abs <- abs(loading_diff)
        int_diff_abs <- abs(intercept_diff)
        maxind <- which(load_diff_abs == max((load_diff_abs)), arr.ind = TRUE)
        maxind2 <- which(int_diff_abs == max((int_diff_abs)), arr.ind = TRUE)

        if (display == TRUE) {
          cat("Factor loadings for each group ---------------------------------------------------------", "\n", "\n")
          print(loadings)
          cat("Difference in loadings between groups --------------------------------------------------", "\n", "\n")
          print(loading_diff)

          print(paste("The largest difference in loading is that of variable",
                      dimnames(loading_diff)[[1]][maxind[1]],
                      "on factor", dimnames(loading_diff)[[2]][maxind[2]], "for group comparison",
                      dimnames(loading_diff)[[3]][maxind[3]]))
          cat("\n")
          cat("Intercepts for each group --------------------------------------------------------------", "\n", "\n")
          print(intercepts)
          cat("\n")
          cat("Difference in intercepts between groups ------------------------------------------------", "\n", "\n")
          print(intercept_diff)
          cat("\n")
          print(paste("The largest difference in intercept is that of variable",
                      dimnames(intercept_diff)[[1]][maxind2[1]],
                      "for group comparison", dimnames(intercept_diff)[[2]][maxind2[2]]))
        }

      } else {
        loading_invariance <- cfa(model, data, std.lv = TRUE,
                                  missing = "fiml", group = group,
                                  group.equal = c("loadings"), estimator = "MLR")

        load_check <- semTools::compareFit(Loadings_Invariant = loading_invariance,
                                           Loadings_Free = cfa)
        summary.FitDiff(load_check)

        cat("Loading Invariance Interpretation -----------------------------------------------------------", "\n", "\n")

        if (pvalue_chisq_diff1 >= 0.05){
          cat("The hypothesis of perfect loading invariance *is not* rejected
           according to the Chi-Square difference test statistics because the p-value
           is larger or equal to 0.05.", "\n", "\n")
        } else {
          cat("The hypothesis of perfect loading invariance *is* rejected according
           to the Chi-Square difference test statistics because the p-value is
           smaller than 0.05", "\n", "\n")
        }

        # CFI  ----------------------------------------------------------------
        if (diff_cfi_scaled1 <= .01) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
           according to the CFI because the difference in CFI value is smaller than
           or equal to 0.01.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
           rejected according to the CFI because the difference in CFI value is
           larger than 0.01.", "\n", "\n")
        }

        # RMSEA ---------------------------------------------------------------
        if(diff_rmsea_scaled1 <= 0.015) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
           according to the RMSEA because the difference in RMSEA value is smaller
           than or equal to 0.015.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
           rejected according to the RMSEA because the difference in RMSEA value is
           larger than 0.015.", "\n", "\n")
        }

        # SRMR ----------------------------------------------------------------
        if (diff_srmr_scaled1 <= .030) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
           according to the SRMR because the difference in SRMR value is smaller than
           or equal to 0.030.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
           rejected according to the SRMR because the difference in SRMR value is
           larger than 0.030.", "\n", "\n")
        }

        t1 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                         as.numeric(levels(pull(data, group))[1]),"`$lambda)")
        )))

        loadings <- array(NA, dim = c(t1[1], t1[2],length(levels(pull(data, group)))))
        group_upd <- as.numeric(pull(data, group))

        for (i in min(group_upd):max(group_upd)){
          eval(parse(text = c(
            paste0("loadings[,,", i, "] <-", "inspect(cfa, \"est\")$`",
                   as.numeric(levels(pull(data, group))[i]), "`$lambda")
          )))
        }

        rowlab <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$lambda)[[1]]")
        )))

        collab <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$lambda)[[2]]")
        )))

        dimnames(loadings) <- list(rowlab, collab,levels(pull(data, group)))

        combis <- as.matrix(combn(unique(group_upd), 2))

        loading_diff <- array(NA, dim = c(t1[1], t1[2], ncol(combis)))
        label_diff <- rep(NA, ncol(combis))

        for (j in 1:ncol(combis)){
          loading_diff[,,j] <- loadings[,,combis[1,j]] - loadings[,,combis[2,j]]
          label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                  levels(pull(data, group))[combis[2,j]])
        }

        dimnames(loading_diff) <- list(rowlab, collab, label_diff)
        load_diff_abs <- abs(loading_diff)
        maxind <- which(load_diff_abs == max((load_diff_abs)), arr.ind = TRUE)

        if (display == TRUE) {
          cat("Factor loadings for each group ---------------------------------------------------------", "\n", "\n")
          print(loadings)
          cat("Difference in loadings between groups --------------------------------------------------", "\n", "\n")
          print(loading_diff)

          print(paste("The highest loading is that of variable", dimnames(loadings)[[1]][maxind[1]],
                      "on factor", dimnames(loadings)[[2]][maxind[2]], "in group",
                      dimnames(loadings)[[3]][maxind[3]]))
        }
        cat("\n")
        cat("Intercept Invariance Interpretation ------------------------------------------------------", "\n", "\n")

        print("Before testing intercept invariance, (partial) loading invariance should hold, which is not the case.")
      }
    }
  }
  if (estimator == "ML") {
    if (intercept == FALSE) {
      if (missing == FALSE) {
        cfa <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                   group = group, estimator = "ML")

        loading_invariance <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                                  group = group,
                                  group.equal = c("loadings"), estimator = "ML")
      } else {
        cfa <- cfa(model, data, std.lv = TRUE, missing  = "fiml",
                   group = group, estimator = "ML")

        loading_invariance <- cfa(model, data, std.lv = TRUE,
                                  missing = "fiml", group = group,
                                  group.equal = c("loadings"), estimator = "ML")
      }

      load_check <- semTools::compareFit(Loadings_Invariant = loading_invariance,
                                         Loadings_Free = cfa)

      summary.FitDiff(load_check)

      pvalue_chisq_diff <- as.numeric(
        anova(loading_invariance, cfa)$"Pr(>Chisq)"[2])

      diff_cfi <- fitMeasures(loading_invariance)[["cfi"]] -
        fitMeasures(cfa)[["cfi"]]
      diff_rmsea <- fitMeasures(loading_invariance)[["rmsea"]] -
        fitMeasures(cfa)[["rmsea"]]
      diff_srmr <- fitMeasures(loading_invariance)[["srmr"]] -
        fitMeasures(cfa)[["srmr"]]

      # Chi-Square ------------------------------------------------------------
      cat("Loading Invariance Interpretations---------------------------------------------------------", "\n", "\n")

      if (pvalue_chisq_diff >= 0.05) {
        cat("The hypothesis of perfect loading invariance *is not* rejected
         according to the Chi-Square difference test statistics because the p-value
         is larger or equal to 0.05.", "\n", "\n")
      } else {
        cat("The hypothesis of perfect loading invariance *is* rejected according
         to the Chi-Square difference test statistics because the p-value is
         smaller than 0.05", "\n", "\n")
      }

      # CFI  ------------------------------------------------------------------
      if (diff_cfi <= .01) {
        cat("The hypothesis of approximate loading invariance *is not* rejected
         according to the CFI because the difference in CFI value is smaller than
         or equal to 0.01.", "\n", "\n")
      } else {
        cat("The hypothesis of approximate approximate loading invariance *is*
         rejected according to the CFI because the difference in CFI value is
         larger than 0.01.", "\n", "\n")
      }

      # RMSEA -----------------------------------------------------------------
      if (diff_rmsea <= 0.015) {
        cat("The hypothesis of approximate loading invariance *is not* rejected
         according to the RMSEA because the difference in RMSEA value is smaller
         than or equal to 0.015.", "\n", "\n")
      } else {
        cat("The hypothesis of approximate approximate loading invariance *is*
         rejected according to the RMSEA because the difference in RMSEA value is
         larger than 0.015.", "\n", "\n")
      }

      # SRMR ------------------------------------------------------------------
      if (diff_srmr < .030) {
        cat("The hypothesis of approximate loading invariance *is not* rejected
         according to the SRMR because the difference in SRMR value is smaller than
         or equal to 0.030.", "\n", "\n")
      } else {
        cat("The hypothesis of approximate approximate loading invariance *is*
         rejected according to the SRMR because the difference in SRMR value is
         larger than 0.030.", "\n", "\n")
      }

      t1 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                       as.numeric(levels(pull(data, group))[1]),"`$lambda)")
      )))

      loadings <- array(NA, dim = c(t1[1], t1[2],length(levels(pull(data, group)))))
      group_upd <- as.numeric(pull(data, group))

      for (i in min(group_upd):max(group_upd)){
        eval(parse(text = c(
          paste0("loadings[,,", i, "] <-", "inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[i]), "`$lambda")
        )))
      }

      rowlab <- eval(parse(text = c(
        paste0("dimnames(inspect(cfa, \"est\")$`",
               as.numeric(levels(pull(data, group))[1]),"`$lambda)[[1]]")
      )))

      collab <- eval(parse(text = c(
        paste0("dimnames(inspect(cfa, \"est\")$`",
               as.numeric(levels(pull(data, group))[1]),"`$lambda)[[2]]")
      )))

      dimnames(loadings) <- list(rowlab, collab, levels(pull(data, group)))

      combis <- as.matrix(combn(unique(group_upd), 2))

      loading_diff <- array(NA, dim = c(t1[1], t1[2], ncol(combis)))
      label_diff <- rep(NA, ncol(combis))

      for (j in 1:ncol(combis)){
        loading_diff[,,j] <- loadings[,,combis[1,j]] - loadings[,,combis[2,j]]
        label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                levels(pull(data, group))[combis[2,j]])
      }

      dimnames(loading_diff) <- list(rowlab, collab, label_diff)
      load_diff_abs <- abs(loading_diff)
      maxind <- which(load_diff_abs == max((load_diff_abs)), arr.ind = TRUE)

      if (display == TRUE) {
        cat("Factor loadings for each group ---------------------------------------------------------", "\n", "\n")
        print(loadings)
        cat("Difference in loadings between groups --------------------------------------------------", "\n", "\n")
        print(loading_diff)

        print(paste("The largest difference in loading is that of variable",
                    dimnames(loading_diff)[[1]][maxind[1]],"on factor",
                    dimnames(loading_diff)[[2]][maxind[2]], "for group comparison",
                    dimnames(loading_diff)[[3]][maxind[3]]))
      }

    } else {
      if (missing == FALSE){
        cfa <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                   group = group, estimator = "ML")

        loading_invariance <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                                  group = group,
                                  group.equal = c("loadings"), estimator = "ML")

        intercept_invariance <- cfa(model, data, std.lv = TRUE, meanstructure = TRUE,
                                    group = group,
                                    group.equal = c("loadings", "intercepts"),
                                    estimator = "ML")

      } else {
        cfa <- cfa(model, data, std.lv = TRUE, missing  = "fiml",
                   group = group, estimator = "ML")

        loading_invariance <- cfa(model, data, std.lv = TRUE,
                                  missing = "fiml", group = group,
                                  group.equal = c("loadings"), estimator = "ML")

        intercept_invariance <- cfa(model, data, std.lv = TRUE,
                                    missing = "fiml", group = group,
                                    group.equal = c("loadings", "intercepts"),
                                    estimator = "ML")
      }


      pvalue_chisq_diff1 <- as.numeric(
        anova(intercept_invariance, loading_invariance, cfa)$"Pr(>Chisq)"[2])

      diff_cfi1 <- fitMeasures(loading_invariance)[["cfi"]] -
        fitMeasures(cfa)[["cfi"]]
      diff_rmsea1 <- fitMeasures(loading_invariance)[["rmsea"]] -
        fitMeasures(cfa)[["rmsea"]]
      diff_srmr1 <- fitMeasures(loading_invariance)[["srmr"]] -
        fitMeasures(cfa)[["srmr"]]

      SC2 <- sum(as.numeric(pvalue_chisq_diff1 >.05),
                 as.numeric(diff_cfi_scaled1 <= .01),
                 as.numeric(diff_rmsea_scaled1 <= 0.015),
                 as.numeric(diff_srmr_scaled1 <= 0.030))

      if (SC2 >= 3){
        int_check <- semTools::compareFit(Intercept_Invariant =
                                            intercept_invariance, Loadings_Invariant = loading_invariance,
                                          Loadings_Intercept_Free = cfa)

        summary.FitDiff(int_check)

        pvalue_chisq_diff2 <- as.numeric(
          anova(intercept_invariance, loading_invariance, cfa)$"Pr(>Chisq)"[3])

        diff_cfi2 <- fitMeasures(intercept_invariance)[["cfi"]] -
          fitMeasures(loading_invariance)[["cfi"]]
        diff_rmsea2 <- fitMeasures(intercept_invariance)[["rmsea"]] -
          fitMeasures(loading_invariance)[["rmsea"]]
        diff_srmr2 <- fitMeasures(intercept_invariance)[["srmr"]] -
          fitMeasures(loading_invariance)[["srmr"]]

        # Chi-Square ----------------------------------------------------------
        cat("Interpretations-------------------------------------------------------
         --------", "\n", "\n")

        cat("Loading Invariance----------------------------------------------------
         --------", "\n", "\n")

        if (pvalue_chisq_diff1 >= 0.05){
          cat("The hypothesis of perfect loading invariance *is not* rejected
           according to the Chi-Square difference test statistics because the p-value
           is larger or equal to 0.05.", "\n", "\n")
        } else {
          cat("The hypothesis of perfect loading invariance *is* rejected according
           to the Chi-Square difference test statistics because the p-value is
           smaller than 0.05", "\n", "\n")
        }

        # CFI  ----------------------------------------------------------------
        if (diff_cfi1 <= .01) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
           according to the CFI because the difference in CFI value is smaller than
           or equal to 0.01.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
           rejected according to the CFI because the difference in CFI value is
           larger than 0.01.", "\n", "\n")
        }

        # RMSEA ---------------------------------------------------------------
        if(diff_rmsea1 <= 0.015) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
           according to the RMSEA because the difference in RMSEA value is smaller
           than or equal to 0.015.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
           rejected according to the RMSEA because the difference in RMSEA value is
           larger than 0.015.", "\n", "\n")
        }

        # SRMR ----------------------------------------------------------------
        if (diff_srmr1 <= .030) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
           according to the SRMR because the difference in SRMR value is smaller than
           or equal to 0.030.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
           rejected according to the SRMR because the difference in SRMR value is
           larger than 0.030.", "\n", "\n")
        }

        cat("Intercept Invariance----------------------------------------------------
         --------", "\n", "\n")

        if (pvalue_chisq_diff2 >= 0.05){
          cat("The hypothesis of perfect intercept invariance *is not* rejected
           according to the Chi-Square difference test statistics because the p-value
           is larger or equal to 0.05.", "\n", "\n")
        } else {
          cat("The hypothesis of perfect intercept invariance *is* rejected
           according to the Chi-Square difference test statistics because the p-value
           is smaller than 0.05", "\n", "\n")
        }

        if (diff_cfi2 <= .01) {
          cat("The hypothesis of approximate intercept invariance *is not*
           rejected according to the CFI because the difference in CFI value is
           smaller than or equal to 0.01.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate intercept invariance *is*
           rejected according to the CFI because the difference in CFI value is
           larger than 0.01.", "\n", "\n")
        }

        if (diff_rmsea2 <= 0.015) {
          cat("The hypothesis of approximate intercept invariance *is not* rejected
           according to the RMSEA because the difference in RMSEA value is smaller
           than or equal to 0.015.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate intercept invariance *is*
           rejected according to the RMSEA because the difference in RMSEA value is
           larger than 0.015.", "\n", "\n")
        }

        if (diff_srmr2 <= .030) {
          cat("The hypothesis of approximate intercept invariance *is not* rejected
           according to the SRMR because the difference in SRMR value is smaller than
           or equal to 0.030.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate intercept invariance *is*
           rejected according to the SRMR because the difference in SRMR value is
           larger than 0.030.", "\n", "\n")
        }

        t1 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                         as.numeric(levels(pull(data, group))[1]),"`$lambda)")
        )))

        loadings <- array(NA, dim = c(t1[1], t1[2],length(levels(pull(data, group)))))
        group_upd <- as.numeric(pull(data, group))

        for (i in min(group_upd):max(group_upd)){
          eval(parse(text = c(
            paste0("loadings[,,", i, "] <-", "inspect(cfa, \"est\")$`",
                   as.numeric(levels(pull(data, group))[i]), "`$lambda")
          )))
        }

        rowlab <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$lambda)[[1]]")
        )))

        collab <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$lambda)[[2]]")
        )))

        dimnames(loadings) <- list(rowlab, collab, levels(pull(data, group)))

        combis <- as.matrix(combn(unique(group_upd), 2))

        loading_diff <- array(NA, dim = c(t1[1], t1[2], ncol(combis)))
        label_diff <- rep(NA, ncol(combis))

        for (j in 1:ncol(combis)){
          loading_diff[,,j] <- loadings[,,combis[1,j]] - loadings[,,combis[2,j]]
          label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                  levels(pull(data, group))[combis[2,j]])
        }

        dimnames(loading_diff) <- list(rowlab, collab, label_diff)

        t2 <- eval(parse(text = c(paste0("dim(inspect(cfa, \"est\")$`",
                                         as.numeric(levels(pull(data, group))[1]),"`$nu)")
        )))

        intercepts <- array(NA, dim = (c(t2[1], length(levels(pull(data, group))))))
        group_upd <- as.numeric(pull(data, group))

        for (i in min(group_upd):max(group_upd)){
          eval(parse(text = c(
            paste0("intercepts[,", i, "] <-", "inspect(cfa, \"est\")$`",
                   as.numeric(levels(pull(data, group))[i]), "`$nu")
          )))
        }

        rowlab2 <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$nu)[[1]]")
        )))

        dimnames(intercepts) <- list(rowlab2, levels(pull(data, group)))

        combis <- as.matrix(combn(unique(group_upd), 2))

        intercept_diff <- array(NA, dim = c(t2[1], ncol(combis)))
        label_diff <- rep(NA, ncol(combis))

        for (j in 1:ncol(combis)){
          intercept_diff[,j] <- intercepts[,combis[1,j]] - intercepts[,combis[2,j]]
          label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                  levels(pull(data, group))[combis[2,j]])
        }

        dimnames(intercept_diff) <- list(rowlab2, label_diff)

        load_diff_abs <- abs(loading_diff)
        int_diff_abs <- abs(intercept_diff)
        maxind <- which(load_diff_abs == max((load_diff_abs)), arr.ind = TRUE)
        maxind2 <- which(int_diff_abs == max((int_diff_abs)), arr.ind = TRUE)

        if (display == TRUE) {
          cat("Factor loadings for each group ---------------------------------------------------------", "\n", "\n")
          print(loadings)
          cat("\n", "\n", "Difference in loadings between groups --------------------------------------", "\n", "\n")
          print(loading_diff)

          print(paste("The largest difference in loading is that of variable",
                      dimnames(loading_diff)[[1]][maxind[1]],
                      "on factor", dimnames(loading_diff)[[2]][maxind[2]], "for group comparison",
                      dimnames(loading_diff)[[3]][maxind[3]]))
          cat("\n")
          cat("Intercepts for each group --------------------------------------------------------------", "\n", "\n")
          print(intercepts)
          cat("\n")
          cat("Difference in intercepts between groups ------------------------------------------------", "\n", "\n")
          print(intercept_diff)
          cat("\n")
          print(paste("The largest difference in intercept is that of variable",
                      dimnames(intercept_diff)[[1]][maxind2[1]],
                      "for group comparison", dimnames(intercept_diff)[[2]][maxind2[2]]))
        }
      } else {
        load_check <- semTools::compareFit(Loadings_Invariant = loading_invariance,
                                           Loadings_Free = cfa)

        summary.FitDiff(load_check)

        cat("Loading Invariance Interpretation ------------------------------------------------------", "\n", "\n")

        if (pvalue_chisq_diff1 >= 0.05){
          cat("The hypothesis of perfect loading invariance *is not* rejected
             according to the Chi-Square difference test statistics because the p-value
             is larger or equal to 0.05.", "\n", "\n")
        } else {
          cat("The hypothesis of perfect loading invariance *is* rejected according
             to the Chi-Square difference test statistics because the p-value is
             smaller than 0.05", "\n", "\n")
        }

        # CFI  --------------------------------------------------------------
        if (diff_cfi1 <= .01) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
             according to the CFI because the difference in CFI value is smaller than
             or equal to 0.01.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
             rejected according to the CFI because the difference in CFI value is
             larger than 0.01.", "\n", "\n")
        }

        # RMSEA -------------------------------------------------------------
        if(diff_rmsea1 <= 0.015) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
             according to the RMSEA because the difference in RMSEA value is smaller
             than or equal to 0.015.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
             rejected according to the RMSEA because the difference in RMSEA value is
             larger than 0.015.", "\n", "\n")
        }

        # SRMR --------------------------------------------------------------
        if (diff_srmr1 <= .030) {
          cat("The hypothesis of approximate loading invariance *is not* rejected
             according to the SRMR because the difference in SRMR value is smaller than
             or equal to 0.030.", "\n", "\n")
        } else {
          cat("The hypothesis of approximate approximate loading invariance *is*
             rejected according to the SRMR because the difference in SRMR value is
             larger than 0.030.", "\n", "\n")
        }

        t1 <- eval(parse(text = c(
          paste0("dim(inspect(cfa, \"est\")$`",
                 as.numeric(levels(data$group)[1]),"`$lambda)")
        )))

        loadings <- array(NA, dim = c(t1[1], t1[2],length(levels(pull(data, group)))))
        group_upd <- as.numeric(pull(data, group))

        for (i in min(group_upd):max(group_upd)){
          eval(parse(text = c(
            paste0("loadings[,,", i, "] <-", "inspect(cfa, \"est\")$`",
                   as.numeric(levels(pull(data, group))[i]), "`$lambda")
          )))
        }

        rowlab <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$lambda)[[1]]")
        )))

        collab <- eval(parse(text = c(
          paste0("dimnames(inspect(cfa, \"est\")$`",
                 as.numeric(levels(pull(data, group))[1]),"`$lambda)[[2]]")
        )))

        dimnames(loadings) <- list(rowlab, collab, levels(pull(data, group)))

        combis <- as.matrix(combn(unique(group_upd), 2))

        loading_diff <- array(NA, dim = c(t1[1], t1[2], ncol(combis)))
        label_diff <- rep(NA, ncol(combis))

        for (j in 1:ncol(combis)){
          loading_diff[,,j] <- loadings[,,combis[1,j]] - loadings[,,combis[2,j]]
          label_diff[j] <- paste0(levels(pull(data, group))[combis[1,j]], "-",
                                  levels(pull(data, group))[combis[2,j]])
        }

        dimnames(loading_diff) <- list(rowlab, collab, label_diff)

        load_diff_abs <- abs(loading_diff)
        maxind <- which(load_diff_abs == max((load_diff_abs)), arr.ind = TRUE)

        if (display == TRUE) {
          cat("Factor loadings for each group ---------------------------------------------------------", "\n", "\n")
          print(loadings)
          cat("Difference in loadings between groups --------------------------------------------------", "\n", "\n")
          print(loading_diff)

          print(paste("The largest difference in loading is that of variable",
                      dimnames(loading_diff)[[1]][maxind[1]],"on factor",
                      dimnames(loading_diff)[[2]][maxind[2]], "for group comparison",
                      dimnames(loading_diff)[[3]][maxind[3]]))
        }

        cat("Intercept Invariance Interpretation --------------------------------------------------------", "\n", "\n")

        print("Before testing intercept invariance, (partial) loading invariance should hold, which is not the case.")
      }
    }
  }
}
