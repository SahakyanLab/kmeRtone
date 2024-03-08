#' Get Virus Metadata Resource (VMR) from International Committee on Taxonomy
#'   of Viruses (ICTV)
#'
#' Always get the latest VMR table, so no argument.
#'
#' @return Virus Metadata Resource data.table.
#'
#' @importFrom data.table setDT
#' @importFrom openxlsx read.xlsx
#' 
#' @export
getICTVvirusMetadataResource <- function() {

  vmr.req <- curl::curl_fetch_memory("https://ictv.global/vmr")

  vmr.urls <- rawToChar(vmr.req$content)
  vmr.urls <- stri_split_fixed(vmr.urls, "\n")[[1]]
  vmr.urls <- vmr.urls[stri_detect_fixed(vmr.urls, "xlsx")]
  vmr.filenames <- stri_extract_first_regex(vmr.urls, "VMR_.+\\.xlsx")
  vmr.urls <- stri_extract_first_regex(vmr.urls, "href.+?>") |>
    stri_extract_first_regex('".+?"') |> stri_replace_all_fixed('"', "")
  vmr.urls <- paste0("https://ictv.global", vmr.urls)

  vmr <- openxlsx::read.xlsx(vmr.urls[1], sheet = 1)
  setDT(vmr)

  return(vmr)
}
