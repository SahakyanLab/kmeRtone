#' A generic function to get Ensembl data persistently from a URL. This is an
#'    internal function used by other getEnsemblXXX functions.
#' 
#' Error is handled based on their rule as set out at
#' https://github.com/Ensembl/ensembl-rest/wiki/HTTP-Response-Codes
#' 
#' @param url Pre-built Ensembl REST URL.
#' @param handle Curl handle object.
#' @param max.attempt Maximum number of attempt. Default to 5.
#' @return A parsed json data. Could be data.frame or list of lists.
#' 
getEnsemblData <- function(url, handle, max.attempt=5) {
  
  attempt.no <- 1
  while (TRUE) {
    
    response <<- curl::curl_fetch_memory(url, handle)
    
    if (response$status_code == 200) {
      
      data <- fromJSON(rawToChar(response$content))
      break
      
    } else if (response$status_code == 400) {
      
      stop("Bad Request (400) ", rawToChar(response$content))
      
    } else if (response$status_code == 403) {
      
      msg <- paste0("403 (Forbidden) You are submitting far too many requests ",
                    "and have been temporarily forbidden access to the ",
                    "service. Wait and retry with a maximum of 15 requests per",
                    " second.")
      stop(msg, " ", rawToChar(response$content))
      
    } else if (response$status_code == 404) {
      
      mesg <- paste0("404 (Not Found) Badly formatted request. Check your URL")
      stop(msg, " ", rawToChar(response$content))
      
    } else if (response$status_code == 408) {
      
      msg <- paste0("408 (Timeout) The request was not processed in time. Wait",
                    " and retry later")
      
      if (attempt.no < max.attempt) Sys.sleep(attempt.no * 5)
      
    } else if (response$status_code == 429) {
      
      headers <- curl::parse_headers_list(response$headers)
      wait.time <- headers$`Retry-After` |> as.numeric()
      msg <- paste0("429 (Too Many Requests) You have been rate-limited; wait ",
                    "and retry. Waiting time is ", wait.time, "seconds.")
      message(msg)
      Sys.sleep(wait.time)

    } else if (response$status_code == 503) {
      
      msg <- paste0("503 (Service Unavailable) The service is temporarily down",
                    "; retry after a pause")
      message(msg)
      if (attempt.no < max.attempt) Sys.sleep(attempt.no * 5)
      
    }
    
    if (attempt.no == max.attempt) {
      
      cat(rawToChar(response$headers), "\n")
      cat(rawToChar(response$content), "\n")
      prompt.msg <- paste0(max.attempt, " attempts have been made. ",
                           "Retry? (y/n/s): ")
      
      while (TRUE) {
        prompt <- readline(prompt.msg)
        if (prompt == "y") {
          Sys.sleep(attempt.no * 5)
          max.attempt <- max.attempt * 2
        } else if (prompt == "n") {
          break
        } else if (prompt == "s") {
          stop("R is stopped as requested.")
        }
      }
    }
    
    attempt.no <- attempt.no + 1
  }
  
  
  return(data)
}