#' Build control regions
#'
#' @param case Case in Coordinate class object format.
#' @param ctrl.rel.pos Control relative position.
#' @param genome Genome class object.
#' @param output.path Output directory path to save control coordinate.
#' @param verbose Boolean. Default is TRUE and will print progress updates. 
#'
#' @return Control in Coordinate class object format.
#'
#' @importFrom data.table fwrite
#' @importFrom future nbrOfWorkers
#' @importFrom progressr progressor
#' 
#' @export
buildControl <- function(case, ctrl.rel.pos, genome, output.path="control/",
                         verbose=TRUE) {

  # Error checking
  if (length(list.files(output.path)) > 0)
    warning("Output directory is not empty.")

  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  ncpu <- future::nbrOfWorkers()

  p <- progressr::progressor(along = case$chr_names)

  future.apply::future_lapply(case$chr_names, function(chr.name) {

    p(chr.name)

    # Print message
    if (verbose) {
      t1 <- Sys.time()
      msg <- paste0("Building control regions of ", chr.name)
      dots <- rep(".", 40 - nchar(msg)) |> paste(collapse = "")
      cat(paste0(msg, dots, if(ncpu > 1) "\n"))
    }

    # Initialize control region
    control <- case[chr.name, state = "case"][, .(
      start = c(start - ctrl.rel.pos[2],
                (if(!is.null(case$single_len)) start + case$single_len - 1 else
                  end) + ctrl.rel.pos[1]),
      end = c(start - ctrl.rel.pos[1],
              (if(!is.null(case$single_len)) (start + case$single_len - 1) else
                end) + ctrl.rel.pos[2])
    )]
    control <- trimCoordinate(control, seq.len = stri_length(genome[chr.name]))
    control <- mergeCoordinate(control)

    # Define buffer region as outermost case regions + ctrl.rel.pos[1]
    buffer <- case[chr.name, state = "case"][
      , .(start = start - ctrl.rel.pos[1],
          end = (if(!is.null(case$single_len)) start + case$single_len - 1 else
            end) + ctrl.rel.pos[1])]
    buffer <- trimCoordinate(buffer, seq.len = stri_length(genome[chr.name]))
    buffer <- mergeCoordinate(buffer)

    # Final control coordinate
    control <- removeRegion(control, region = buffer)

    # Save coordinates
    fwrite(control, paste0(output.path, "/", chr.name, ".csv"),
           showProgress = FALSE)

    # Print message
    if (verbose) {
      t <- Sys.time() - t1
      cat(paste0(if(ncpu > 1) paste0(msg, dots), "DONE! -- ", round(t[[1]], 2),
                 " ", attr(t, "units"), "\n"))
    }

  }, future.seed = NULL)

  # Build Coordinate object
  control <- Coordinate$new(root.path = output.path,
                            is.strand.sensitive = FALSE)

  return(control)
}
