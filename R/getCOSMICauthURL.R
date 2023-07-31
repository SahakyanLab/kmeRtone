#' Get COSMIC authenticated URL.
#' 
#' To access the data for non-commercial usage, you must register with the
#'    COSMIC. This function fetch the authenticated URL from the public URL
#'    given by the COSMIC website.
#'
#' @param email Email registered with COSMIC.
#' @param password Password.
#' @param url URL given by COSMIC website to access data.
#'
#' @return An authenticated URL, valid for 1 hour access.
#'
#' @export
getCOSMICauthURL <- function(email, password, url) {
  
  id <- paste0(email, ":", password)
  raw.id <- charToRaw(id)
  auth.header <- paste("Basic", jsonlite::base64_enc(raw.id))
  
  h <- curl::new_handle()
  curl::handle_setheaders(h, "Authorization" = auth.header)
  
  response <- curl::curl_fetch_memory(url, handle = h)
  auth.url <- jsonlite::fromJSON(rawToChar(response$content))$url
  
  return(auth.url)
}