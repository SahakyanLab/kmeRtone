.onLoad <- function(libname, pkgname) {
  if (Sys.getenv("R_KMERTONE_DATA") == "")
    Sys.setenv(R_KMERTONE_DATA = "~/kmertone_data")
  assign("kmertone.data.path", Sys.getenv("R_KMERTONE_DATA"), envir = topenv())
}
