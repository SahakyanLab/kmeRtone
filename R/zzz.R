.onLoad <- function(libname, pkgname) {
  if (Sys.getenv("R_kmeRtone_DATA") == "")
    Sys.setenv(R_kmeRtone_DATA = "~/kmeRtone_data")
  assign("kmeRtone.data.path", Sys.getenv("R_kmeRtone_DATA"), envir = topenv())
}
