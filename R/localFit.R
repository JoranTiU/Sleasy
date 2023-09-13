localFit <- function(x){

  local_misfit <- abs(inspect(x, "resid")$cov)
  max_local_misfit <- which(local_misfit == max(local_misfit), arr.ind = TRUE)[1,]

  return(list(local_misfit = local_misfit,
              max_misfit = local_misfit[as.numeric(max_local_misfit[1]),as.numeric(max_local_misfit[2])]))

}
