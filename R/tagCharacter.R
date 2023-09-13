tagCharacter <- function(vec, char = "\u2020", minvalue = NA,
                         print_integer = FALSE, nd = 3L) {
  char <- if (is.null(char)) as.character(NA) else as.character(char)[1]
  if (nchar(char) != 1L) {
    message('Only a single character can be used to tag= preferred models, so ',
            'no tags were added.  To omit tags, specify tag=NULL or tag=NA.')
    char <- as.character(NA)
  }
  if (print_integer) {
    vec <- noLeadingZero(vec, fmt = "%.0f", nd = nd)
  } else if (is.na(minvalue)) {
    vec <- noLeadingZero(vec, fmt = paste0("%.", nd, "f"), nd = nd)
  } else {
    target <- if (minvalue) min(vec, na.rm = TRUE) else max(vec, na.rm = TRUE)
    tag <- rep(" ", length(vec))
    if (!is.na(char)) tag[vec == target] <- char
    vec <- noLeadingZero(vec, fmt = paste0("%.", nd, "f"), nd = nd)
    vec <- paste0(vec, tag)
  }

  vec
}
