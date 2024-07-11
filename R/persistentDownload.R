#' Download file until successful
#'
#' If download failed, it will be repeated until max attempt reached.
#'
#' @param file.url File uniform resource locator.
#' @param output.name Output name.
#' @param max.attempt Maximum number of attempt. Default is 5.
#' @param user.invoke If number of attempt is reached, ask user whether to keep
#'    trying. Default is TRUE to invoke the prompt.
#' @param header A named list or vector of curl header.
#'
#' @return A downloaded file.
#' 
#' @importFrom curl new_handle handle_setheaders curl_fetch_disk
#'
#' @export
persistentDownload <- function(file.url, output.name, max.attempt = 5,
                               user.invoke = TRUE, header) {

  dir.create(dirname(output.name), showWarnings = FALSE, recursive = TRUE)
  
  h <- curl::new_handle()
  if (!missing(header)) curl::handle_setheaders(h, .list = header)
  
  attempt <- 1
  while (TRUE) {

    # If error, response is a list with 2 named elements: message and call
    response <- tryCatch(curl::curl_fetch_disk(file.url, output.name, h),
                         error = function(e) e)
    
    if (sum(response$status_code) == "200") {

      break

    } else if (attempt <= max.attempt) {

      Sys.sleep(attempt * 3)

    } else {

      if (!is.null(response$message)) message(response$message)
      message(paste("Download status: ",
          if (!is.null(response$status_code)) response$status_code
          else "unknown",
          "\n",
          "Current time: ", paste(Sys.time()), "\n",
          "File URL: ", file.url, "\n",
          "Output name: ", output.name, "\n",
          "Download failed after ", attempt-1, " attempts.\n", sep = ""))
      user.response <- NA
      while(!user.response %in% c("y", "n", "s")) {
        user.response <- readline("Keep trying? (y/n/s): ")
      }
      if (user.response == "n") {
        unlink(output.name)
        break
      } else if (user.response == "y") {
        max.attempt <- max.attempt * 2
      } else if (user.response == "s") {
        unlink(output.name)
        stop("Stop R as requested.")
      }
    }
    attempt <- attempt + 1
  }

  return(invisible(NULL))
}
