buildRESTurl <- function(url, .list=list(), ...) {
  
  query <- append(.list, list(...))
  query <- paste(names(query), query, sep = "=", collapse = ";")
  rest.url <- paste0(url, "?", query)
  
  return(rest.url)
}