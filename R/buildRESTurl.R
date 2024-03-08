#' Function constructs a URL for a REST API call by appending query parameters.
#' 
#' @param url Base URL of the REST API.
#' @param .list A list of named query parameters.
#' @param ... additional optional arguments
#' @return string of the full REST API URL.
#'
#' @export
buildRESTurl <- function(url, .list=list(), ...) {
  
  query <- append(.list, list(...))
  query <- paste(names(query), query, sep = "=", collapse = ";")
  rest.url <- paste0(url, "?", query)
  
  return(rest.url)
}