#' Get features of given variant IDs.
#'
#' Function retrieves features for given variant IDs from the Ensembl
#' database. It handles requests asynchronously in batches due to server limits
#' and includes options to fetch additional variant information. Error handling
#' for different HTTP response statuses is implemented to manage request failures.
#'
#' @param species Species name or alias (e.g., homo_sapiens, human).
#' @param variant.ids A vector of variant IDs (e.g., rs56116432, COSM476).
#' @param include.genotypes Include genotypes in the response? Default FALSE.
#' @param include.phenotypes Include phenotypes in the response? Default FALSE.
#' @param include.allele.frequencies Include allele frequencies? Default FALSE.
#' @param include.genotype.frequencies Include genotype frequencies? Default FALSE.
#' @param curl.max.con Maximum number of concurrent connections for curl requests.
#'        Default is 100.
#' @param verbose Verbosity level: 1 for error only, 2 for all requests. Default 1.
#' 
#' @return A variant-named list containing lists of variant features.
#'
#' @importFrom curl new_handle handle_setheaders handle_setopt new_pool multi_add 
#'                 multi_run handle_setopt parse_headers_list
#' @importFrom RcppSimdJson fparse
#' 
#' @export
getEnsemblVariantFeatures <- function(species, variant.ids,
                                      include.genotypes=FALSE,
                                      include.phenotypes=FALSE,
                                      include.allele.frequencies=FALSE,
                                      include.genotype.frequencies=FALSE,
                                      curl.max.con=100, verbose=1) {

  # Build REST URL
  url <- paste0("https://rest.ensembl.org/variation/", species)
  includes <- c(include.genotypes, include.phenotypes,
                include.allele.frequencies, include.genotype.frequencies) |>
    as.integer() |> as.list()
  names(includes) <- c("genotypes", "phenotypes", "pops",
                       "population_genotypes")
  url <- buildRESTurl(url, .list = includes)

  # Count error codes in a error-named vector.
  error.codes <- numeric()

  # Wait time for error 429
  wait.429 <- 0

  # Final content for error 400
  content.400 <- character()

  # Final content for error 429
  content.429 <- character()

  # Vector of failed download json
  failed.variant.ids.jsons <- character()

  # Save content result to temporary location.
  contents.dat <- tempfile(fileext = ".dat")
  contents.dat.con <- file(contents.dat, "ab")
  writeBin(charToRaw("{"), contents.dat.con)

  generateCallbackDone <- function(ids.json) {

    callback <- function(response) {

      if (response$status_code == 200) {

        # Replace } bracket with a comma.
        response$content[length(response$content)] <- charToRaw(",")

        # Remove { bracket and save
        writeBin(response$content[-1], contents.dat.con, useBytes = TRUE)

      } else {

        failed.variant.ids.jsons <<- c(failed.variant.ids.jsons,
                                       ids.json)

        error.codes[as.character(response$status_code)] <<-
          sum(1, error.codes[as.character(response$status_code)], na.rm = TRUE)

        if (response$status_code == 429) {
          headers <- curl::parse_headers_list(response$headers)
          wait.429 <<- headers$`Retry-After` |> as.numeric()
        } else if (response$status_code == 400) {
          content.400 <<- rawToChar(response$content)
        } else if (response$status_code == 429) {
          content.429 <<- rawToChar(response$content)
        }

      }
    }
  }

  # Invent error code 999 for callback fail
  cb.fail.code <- "999"
  cb.fail.msg <- NULL
  generateCallbackFail <- function(ids.json) {
    callback <- function(msg) {
      error.codes[cb.fail.code] <<- sum(1, error.codes[cb.fail.code],
                                         na.rm = TRUE)
      cb.fail.msg <<- msg
      failed.variant.ids.jsons <<- c(failed.variant.ids.jsons, ids.json)
    }
    return(callback)
  }

  # Server limit 15 requests per second. At the moment curl has no option for
  # async option to limit the rate; only max concurrent connection.
  # total_con (max total concurrent connections) is default to 100.
  # host_con (max concurrent connections per host) is default to 6.
  # Use "lsof -i" command to monitor establish TCP port.
  pool <- curl::new_pool(host_con = curl.max.con, total_con = curl.max.con)

  # Server limit to 200 variant IDs per request.
  max.id.req <- 200
  batch.size <- ceiling(length(variant.ids) / max.id.req)
  req.num <- 0
  out <- list(done = 0, fail = 0, pending = 0)
  while (TRUE) {

    # Priotize initial requests.
    if (req.num < batch.size) {
      i <- 1 + (req.num * max.id.req)
      variant.ids.json <- paste0('{"ids" : ["',
                                 paste(variant.ids[i:min(length(variant.ids),
                                                         i + max.id.req - 1)],
                                       collapse = '", "'),
                                 '" ]}')

      req.num <- req.num + 1

    # Then additional requests for failed download.
    } else if (length(failed.variant.ids.jsons) > 0) {
      variant.ids.json <- failed.variant.ids.jsons[1]
      failed.variant.ids.jsons <- failed.variant.ids.jsons[-1]
      req.num <- req.num + 1
    } else {
      variant.ids.json <- character()
    }

    # Free up connection to add another one.
    if (req.num >= curl.max.con && (out$pending + 15) >= curl.max.con) {
      n <- ifelse(out$pending >= curl.max.con, out$pending - curl.max.con + 15,
                  curl.max.con - out$pending)
      out <- curl::multi_run(pool = pool, poll = n)
    }

    if (length(variant.ids.json) == 1) {
      if (verbose == 2) {
        cat("Current request number:", req.num, "out of", batch.size, "\n")
      }
      curl::multi_add(
        handle = curl::new_handle(url = url) |>
          curl::handle_setheaders("content-type" = "application/json",
                                  accept = "application/json") |>
          curl::handle_setopt(customrequest = "POST",
                              postfields = variant.ids.json),
        done = generateCallbackDone(variant.ids.json),
        fail = generateCallbackFail(variant.ids.json),
        pool = pool
      )
    }

    # Sleep for 1 second for every 15 request.
    if (req.num %% 15 == 0) {
      out <- curl::multi_run(pool = pool, poll = TRUE)
      Sys.sleep(1)
    }

    # This is to trigger multi_run() when final batch is less than 15 and when
    # repeating download.
    if (req.num == batch.size | "999" %in% names(error.codes) |
         (req.num > batch.size & length(failed.variant.ids.jsons) == 0)) {
      out <- curl::multi_run(pool = pool, poll = TRUE)
    }

    if (length(error.codes) > 0) {

      # Just to see how many errors.
      if (verbose == 1) {
        cat("Current request number:", req.num, "out of", batch.size, "\n")
      }
      if (verbose %in% 1:2) {
        cat("Total errors: ", sum(error.codes), " (",
            round(sum(error.codes) / ceiling(length(variant.ids) / 200) * 100,
                  digits = 2),
            "%)\n", sep = "")
        cat(paste0(names(error.codes), ": ", error.codes, collapse = "\n"),
            "\n")
        if ("400" %in% names(error.codes)) cat("Error 400 content\n",
                                               content.400)
        if ("429" %in% names(error.codes)) cat("Error 429 content\n",
                                               content.429)
        if ("999" %in% names(error.codes)) cat("Error 999 content\n",
                                               cb.fail.msg, "\n")
      }

      if ("429" %in% names(error.codes) &
           (length(wait.429) > 0 && wait.429 > 0)) {
        # wait.429 is empty numeric when request too many at once.
        # Just wait for awhile.
        msg <- paste0("Server max requests. You have been rate-limited; wait ",
                      "and retry. Waiting time is ", wait.429, "seconds.")

        # If I got error here, wait.429 is NA which means the server returns a
        # character which cannot be coerce to numeric.
        message(msg)
        Sys.sleep(wait.429)

      } else if (any(c("503", "999") %in% names(error.codes))) {
        Sys.sleep(9)
      } else if ("500" %in% names(error.codes)) {
        Sys.sleep(15)
      } else {
        Sys.sleep(2)
      }

      # Reset error codes
      error.codes <- numeric()

    } else if (req.num >= batch.size && (out$pending == 0 & out$error == 0 &
                 length(failed.variant.ids.jsons) == 0)) {
      break
    }

  }

  # Close appended json content
  close(contents.dat.con)

  # Read and parse json
  contents <- readBin(contents.dat, what = "raw",
                      n = file.info(contents.dat)$size)
  contents[length(contents)] <- charToRaw("}") # replace comma with bracket.

  feature.lists <- RcppSimdJson::fparse(contents)
  #feature.lists <- jsonlite::fromJSON(rawToChar(response$content))

  missing.ids <- variant.ids[!variant.ids %in% names(feature.lists)]
  if (length(missing.ids) > 0) {
    message(length(missing.ids), " variant IDs are failed to be retrieved.")
  }

  return(feature.lists)
}
