keyword2regex <- function(general.keyword, exact.keyword) {
  # Convert keyword to regex expression
  # The string that regex is being applied to is separated by comma
  # e.g. "melanoma, sarcoma, LSD, neck sarcoma"

  # Resolve exact-keyword regex expression
  begin.rgx <- "((^|, )"
  end.rgx <- "($|,))"
  exact.keyword <- CJ(begin.rgx, exact.keyword, end.rgx)
  exact.keyword <- exact.keyword[, do.call(paste0, .SD)]

  # Combine exact and general keywords
  regex.expr <- paste(c(exact.keyword, general.keyword),
                      collapse = "|")

  return(regex.expr)
}
