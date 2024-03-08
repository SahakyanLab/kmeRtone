#' A system2 wrapper. If anything happen, just give me error!
#' 
#' Turn warning to error.
#' 
#' @param command the system command to be invoked, as a character string.
#' @param args a character vector of arguments to `command`.
#' @param stdout,stderr where output to ‘stdout’ or ‘stderr’ should be sent.
#'    Possible values are "", to the R console (the default), `NULL` or `FALSE`
#'    (discard output), `TRUE` (capture the output in a character vector) or a
#'    character string naming a file.
#' @param stdin should input be diverted? "" means the default, alternatively a
#'    character string naming a file. Ignored if input is supplied.
#' @param input if a character vector is supplied, this is copied one string per
#'    line to a temporary file, and the standard input of command is redirected
#'    to the file.
#' @param env character vector of name=value strings to set environment
#'    variables.
#' @param wait a logical (not `NA`) indicating whether the `R` interpreter
#'    should wait for the command to finish, or run it asynchronously. This
#'    will be ignored (and the interpreter will always wait) if stdout = `TRUE`
#'    or stderr = `TRUE.` When running the command asynchronously, no output
#'    will be displayed on the Rgui console in Windows (it will be dropped,
#'    instead).
#' @param minimized,invisible arguments that are accepted on Windows but ignored
#'    on this platform, with a warning.
#' @param timeout timeout in seconds, ignored if 0. This is a limit for the
#'    elapsed time running command in a separate process. Fractions of seconds
#'    are ignored.
#' 
system3 <- function(command, args = character(), stdout = "", stderr = "", 
                    stdin = "", input = NULL, env = character(), wait = TRUE, 
                    minimized, invisible, timeout = 0) {
  
  tryCatch(expr = system2(command, args, stdout, stderr, stdin, input, env,
                          wait, minimized, invisible, timeout),
           error = function(e) stop(e),
           warning = function(w) stop(w))
  
}