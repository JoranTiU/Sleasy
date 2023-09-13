noLeadingZero <- function(vec, fmt, nd = 3L) {
  out <- sprintf(fmt, vec)
  upper.limit <- paste0(".", paste(rep(9, nd - 1L), collapse = ""), "5")
  used <- vec < as.numeric(upper.limit) & vec >= 0
  used[is.na(used)] <- FALSE
  out[used] <- substring(out[used], 2)
  out
}
