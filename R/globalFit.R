globalFit <- function(x){

  chisq_scaled <- NA
  chisq_scaled <- fitMeasures(x)[["chisq.scaled"]]

  if(!is.na(chisq_scaled)){
    df_chisq_scaled <-fitMeasures(x)[["df.scaled"]]
    pvalue_chisq_scaled <- fitMeasures(x)[["pvalue.scaled"]]
    cfi_scaled <- fitMeasures(x)[["cfi.scaled"]]
    rmsea_scaled <- fitMeasures(x)[["rmsea.scaled"]]
    rmsea_lower_scaled <- fitMeasures(x)[["rmsea.ci.lower.scaled"]]
    rmsea_upper_scaled <- fitMeasures(x)[["rmsea.ci.upper.scaled"]]
    srmr_scaled <- fitMeasures(x)[["srmr_mplus"]]
  }

  chisq <- fitMeasures(x)[["chisq"]]
  df_chisq <- fitMeasures(x)[["df"]]
  pvalue_chisq <- fitMeasures(x)[["pvalue"]]
  cfi <- fitMeasures(x)[["cfi"]]
  rmsea <- fitMeasures(x)[["rmsea"]]
  rmsea_lower <- fitMeasures(x)[["rmsea.ci.lower"]]
  rmsea_upper <- fitMeasures(x)[["rmsea.ci.upper"]]
  srmr <- fitMeasures(x)[["srmr"]]
  null_rmsea <- as.numeric(nullRMSEA(x, silent = TRUE))
  #null_rmsea <- 0.100 # For testing purposes



  # normal =====================================================================

  cat("Results------------------------------------------------------------------------",
      "\n","\n")

  if(is.na(chisq_scaled)){
    cat("Chi-Square (" ,df_chisq ,") = ",chisq, " with p-value = ",
        pvalue_chisq, "\n","\n",
        "CFI = ", cfi, "\n","\n",
        "RMSEA = ", rmsea, "; lower bound = " ,rmsea_lower, "; upper bound = ",
        rmsea_upper,"\n","\n",
        "SRMR = ", srmr, "\n","\n",
        sep = "")

    #Chi-Square-----------------------------------------------------------------
    cat("Interpretations---------------------------------------------------------------",
        "\n","\n")

    if(pvalue_chisq >= 0.05){
      cat("The hypothesis of perfect fit *is not* rejected according to the
          Chi-Square test statistics because the p-value is larger or equal
          to 0.05.","\n","\n")
    }else{
      cat("The hypothesis of perfect fit *is* rejected according to the Chi-
          Square test statistics because the p-value is smaller than 0.05",
          "\n","\n")
    }

    #CFI  ----------------------------------------------------------------------
    if(null_rmsea < 0.158){
      cat("The RMSEA of the null-model is less than .158. Therefore the CFI is
          likely not informative as a measure of fit for this model, and will
          not be reported","\n","\n")
    }else{
      if(cfi > 0.9){
        cat("The hypothesis of approximate model fit *is not* rejected according
          to the CFI because the value is larger than 0.9.","\n","\n")
      }else{
        cat("The hypothesis of approximate model fit *is* rejected according to
          the CFI because the value is smaller or equal than 0.9.","\n","\n")
      }
    }

    #RMSEA----------------------------------------------------------------------
    if(rmsea < 0.05){
      if(rmsea_upper < 0.05){
        cat("The hypothesis of approximate model fit *is not* rejected according
         to the RMSEA because the point estimate and the upper-limit of it's
               95% CI are smaller than 0.05.","\n","\n")
      }else{cat("The RMSEA is smaller than 0.05, but the upper-limit of it's
               95% CI is not. Therefore the hypothesis of approximate model
               fit *might be* rejected according to the RMSEA","\n","\n")}
    }else{
      cat("The hypothesis of approximate model fit *is* rejected according
         to the RMSEA because the point estimate is larger or equal to
         0.05.","\n","\n")
    }

    #SRMR-----------------------------------------------------------------------
    if(srmr < 0.08){
      cat("The hypothesis of approximate model fit *is not* rejected according
          to the SRMR because the value is smaller than 0.08.",
          "\n","\n")
    }else{
      cat("The hypothesis of approximate model fit *is* rejected according
          to the SRMR because the value is larger or equal to
          0.08.","\n","\n")
    }
    # non-normal==================================================================
  }else{

    cat("Chi-Square (" ,df_chisq_scaled ,") = ",chisq_scaled, " with p-value
          = ", pvalue_chisq_scaled, "\n","\n",
        "CFI = ", cfi_scaled, "\n","\n",
        "RMSEA = ", rmsea_scaled, "; lower bound = " ,rmsea_lower_scaled, ";
      upper bound = ", rmsea_upper_scaled,"\n","\n",
        "SRMR = ", srmr_scaled, "\n","\n",
        sep = "")

    cat("Interpretations---------------------------------------------------------------",
        "\n","\n")
    #Chi-Square-----------------------------------------------------------------
    if(pvalue_chisq_scaled >= 0.05){
      cat("The hypothesis of perfect fit *is not* rejected according to the
          Chi-Square test statistics because the p-value is larger or equal
          to 0.05.","\n","\n")
    }else{
      cat("The hypothesis of perfect fit *is* rejected according to the Chi-
          Square test statistics because the p-value is smaller than 0.05",
          "\n","\n")
    }

    #CFI  ----------------------------------------------------------------------
    if(null_rmsea < 0.158){
      cat("The RMSEA of the null-model is less than .158. Therefore the CFI is
          likely not informative as a measure of fit for this model, and will
          not be reported","\n","\n")
    }else{
      if(cfi_scaled > 0.9){
        cat("The hypothesis of approximate model fit *is not* rejected according
          to the CFI because the value is larger than 0.9.","\n","\n")
      }else{
        cat("The hypothesis of approximate model fit *is* rejected according to
          the CFI because the value is smaller or equal than 0.9.","\n","\n")
      }
    }

    #RMSEA----------------------------------------------------------------------
    if(rmsea_scaled < 0.05){
      if(rmsea_upper_scaled < 0.05){
        cat("The hypothesis of approximate model fit *is not* rejected according
         to the RMSEA because the point estimate and the upper-limit of it's
               95% CI are smaller than 0.05.","\n","\n")
      }else{cat("The RMSEA is smaller than 0.05, but the upper-limit of it's
               95% CI is not. Therefore the hypothesis of approximate model
               fit *might be* rejected according to the RMSEA","\n","\n")}
    }else{
      cat("The hypothesis of approximate model fit *is* rejected according
         to the RMSEA because the point estimate is larger or equal to
         0.05.","\n","\n")
    }
    #SRMR-----------------------------------------------------------------------
    if(srmr_scaled < 0.08){
      cat("The hypothesis of approximate model fit *is not* rejected according
         to the SRMR because the value is smaller than 0.08.",
          "\n","\n")
    }else{
      cat("The hypothesis of approximate model fit *is* rejected according
         to the SRMR because the value is larger or equal to 0.08.",
          "\n","\n")
    }

  }
}
