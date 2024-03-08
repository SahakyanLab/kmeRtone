#' Get COSMIC authenticated URL.
#' 
#' To access the data for non-commercial usage, you must register with the
#'    COSMIC. This function fetch the authenticated URL from the public URL
#'    given by the COSMIC website.
#'
#' @param email Email registered with COSMIC.
#' @param password Password associated with the registered email.
#' @param url Public URL provided by the COSMIC website for data access.
#'
#' @return Authenticated URL valid for 1-hour access to COSMIC data.
#'
#' @importFrom jsonlite fromJSON base64_enc
#' @importFrom curl new_handle handle_setheaders curl_fetch_memory
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