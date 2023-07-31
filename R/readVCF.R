#' Read VCF file using tabix.
#'
#' Require tabix in PATH
#' VCF manual is referred from https://samtools.github.io/hts-specs/VCFv4.3.pdf
#'
#' @param vcf.file A path to a local or remote tabix-indexed VCF file.
#' @param chr.names Chromosome names.
#' @param starts Start positions.
#' @param ends End positions.
#' @param INFO.filter Parse only filtered INFO ID. Default is to parse all IDs.
#'
#' @return A data.table of VCF.
#'
#' @export
readVCF <- function(vcf.file, chr.names, starts, ends, INFO.filter=NULL) {

  metainfo <- getVCFmetainfo(vcf.file)

  header <- (metainfo[length(metainfo)] |> stri_replace_first_fixed("#", "") |>
               stri_split_fixed("\t"))[[1]]

  info.ids <- metainfo[grepl("##INFO", metainfo)] |>
    stri_extract_first_regex("ID=.+?,") |> stri_replace_all_regex("^ID=|,$", "")

  # Filter INFO field to parse
  if (!is.null(INFO.filter)) {
    missing.ids <- INFO.filter[!INFO.filter %in% info.ids]
    if (length(missing.ids) > 0) {
      stop("These INFO IDs are not available:", missing.ids)
    } else {
      info.ids <- info.ids[info.ids %in% INFO.filter]
    }

    # Update metainfo to only exclude other than filtered INFO
    idx.info <- grep("##INFO", metainfo)
    rgx.info <- paste(INFO.filter, collapse = "|")
    idx.info.filter <- grep(paste0("##INFO.+?ID=(", rgx.info, "),"), metainfo)
    idx.exclude <- idx.info[!idx.info %in% idx.info.filter]
    metainfo <- metainfo[-idx.exclude]
  }

  info.types <- metainfo[grepl("##INFO", metainfo)] |>
    stri_extract_first_regex("Type=.+?,") |>
    stri_replace_all_regex("^.+=|,$", "")

  info.numbers <- metainfo[grepl("##INFO", metainfo)] |>
    stri_extract_first_regex("Number=.+?,") |>
    stri_replace_all_regex("^.+=|,$", "")

  # Read the VCF main content
  if (!missing(chr.names)) {
    tmp.vcf <- tempfile(fileext = ".vcf")
    query.regions <- paste0(chr.names, ":", starts, "-", ends)
    system3("tabix", args = c("-D", vcf.file, query.regions, ">", tmp.vcf))
    if (file.size(tmp.vcf) > 0) {
      vcf.dt <- fread(tmp.vcf, sep = "\t", col.names = header,
                      showProgress = FALSE)
    } else {
      # construct empty table and return

      header <- header[header != "INFO"]

      info.cols <- lapply(seq_along(info.ids), function(i) {

        # Separate to create list column
        if (grepl("^[1-9ARG.]", info.numbers[i])) {
          return(list())
        } else if (info.types[i] %in% c("Integer", "Float")) {
          return(numeric())
        } else if (info.types[i] == "Flag") {
          return(logical())
        }
      })
      names(info.cols) <- paste0("INFO.", info.ids)

      cols <- rep(list(character()), length(header))
      names(cols) <- header
      cols[["POS"]] <- numeric()
      cols <- cols[names(cols) != "INFO"]
      cols <- c(cols, info.cols)
      return(as.data.table(cols))
    }

  } else {
    vcf.dt <- fread(vcf.file, sep = "\t", col.names = header,
                    skip = length(metainfo), showProgress = FALSE)
  }

  # Assign list column for ALT
  vcf.dt[, ALT := stri_split_fixed(ALT, ",")]

  # Expand INFO data to INFO.??? columns
  for (i in seq_along(info.ids)) {

    # Check if there is no more allocation for data.table column
    if (length(vcf.dt) >= truelength(vcf.dt))
      setalloccol(vcf.dt, truelength(vcf.dt) + 1)

    idx <- stri_locate_first_regex(vcf.dt$INFO, paste0("(^|;)", info.ids[i],
                                                       "((=.+?(;|$))|(;|$))"))

    info <- stri_sub(vcf.dt$INFO, idx) |>
      stri_replace_first_regex(".+?=", "") |> stri_replace_last_fixed(";", "")

    # To speed up regex match, remove previously matched regex.
    set(vcf.dt, j = "INFO",
        value = stri_sub_replace(vcf.dt$INFO, idx, replacement = ";",
                                 omit_na = TRUE))

    # Separate to create list column
    if (grepl("^[1-9ARG.]", info.numbers[i])) {
      # "." is empty field. Change to empty string to coerce to NA w/o warning
      info <- stri_split_fixed(info, ",") |>
        lapply(function(info) {info[info == "."] <- ""; info})
    }

    # Assign numeric
    if (info.types[i] %in% c("Integer", "Float")) {
      if (grepl("^[1-9ARG.]", info.numbers[i])) {
        info <- lapply(info, as.numeric)
      } else {
        info <- as.numeric(info)
      }
    } else if (info.types[i] == "Flag") {
      info <- !is.na(info)
    }

    set(vcf.dt, j = paste0("INFO.", info.ids[i]), value = info)

  }

  set(vcf.dt, j = "INFO", value = NULL)

  # This is enough for my purpose. R6 class would be better.
  class(vcf.dt) <- c("VCF", class(vcf.dt))
  setattr(vcf.dt, "metainfo", metainfo)

  return(vcf.dt)
}