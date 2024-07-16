#' Extract k-mers from a given Coordinate object and Genome objects
#'
#' A k-mer table is initialized and updated in every chromosome-loop operation.
#'      There are 3 modes of extraction. (1) When k is smaller than 9 or k is
#'      larger than 15, the k-mer is extracted in a standard way. A k-mer table
#'      with every possible k-mers is created and updated. (2) For k between
#'      9 and 13, the k-mer sequence is split to half to reduce memory usage
#'      significantly. e.g. ACGTACGTA will become ACGT ACGTA. (3) When k is
#'      larger than 14, k-mers are extracted the same way as (1) but the k-mer
#'      table is grown or expanded for every new k-mer found.
#'
#' @param coor Coordinate class object.
#' @param genome Genome class object.
#' @param k Length of k-mer.
#' @param central.pattern Central pattern of the k-mer, if applicable.
#' @param rm.overlap.region Boolean indicating if overlapping regions should be 
#'    removed. Default is TRUE.
#' @param verbose Boolean indicating if verbose output is enabled.
#'
#' @return A k-mer table with counts for each k-mer.
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom data.table setkey
#' @importFrom progressr progressor
#' @importFrom stringi stri_sub
#' 
#' @export
extractKmers <- function(coor, genome, k, central.pattern=NULL,
                         rm.overlap.region=TRUE, verbose=TRUE) {

  # Error checking
  stopifnot("Coordinate" %in% class(coor))
  stopifnot(any(class(genome) %in% c("NCBI_Genome", "UCSC_Genome")))

  do.standard.extraction <- k <= 8 | k >= 14
  do.split.kmer.extraction <- k >= 9 & k <= 13
  do.growing.kmer.table <- k >= 14

  # Assume small single-length regions to be the case to expand to k-mer.
  expand.to.kmer <- !is.null(coor$single_len) && coor$single_len < k

  # Repetitive function 1
  # Initialize k-mer table
  initialize.kmer.table <- function() {

    if ((do.standard.extraction | do.split.kmer.extraction) &
        !do.growing.kmer.table) {

      kmer.table <- initKmerTable(
        k = k,
        central.pattern = c(central.pattern,
                            if(!is.null(central.pattern) &
                               !coor$is_strand_sensitive)
                              reverseComplement(central.pattern)),
        split.kmer = ifelse(do.split.kmer.extraction, TRUE, FALSE))

    } else if (do.growing.kmer.table) {

      # Empty table
      kmer.table <- initKmerTable(
        split.kmer = ifelse(do.split.kmer.extraction, TRUE, FALSE))

    }

    return(kmer.table)
  }

  # Repetitive function 2
  # A repetitive helper function with exact same condition but with
  #     varied start, end and k because of splitting k-mer counting
  #     method.
  map.kmers <- function(chr.name, start, k, end, rm.trunc.kmer) {
    kmers <- mapKmers(
      seq = genome[chr.name],
      start = start,
      end = end,
      len = if(!expand.to.kmer & !is.null(coor$single_len))
        coor$single_len,
      k = k,
      rm.trunc.kmer = rm.trunc.kmer)
  }

  # MAIN -----------------------------------------------------------------------

  ncpu <- future::nbrOfWorkers()
  kmer.table <- initialize.kmer.table()

  # Environment for updating table row. This is needed because operation "by"
  # in data.table create a new environment for every loop.
  if (do.growing.kmer.table) {
    env <- environment()
  }

  p <- progressr::progressor(along = coor$chr_names)

  kmer.counts <- future.apply::future_lapply(coor$chr_names, function(chr.name){

    p(chr.name)

    # Initialize new k-mer table for every CPU core environment
    if (ncpu > 1) {
      kmer.table <- initialize.kmer.table()

      # Environment for updating row.
      if (do.growing.kmer.table) {
        env <- environment()
      }
    }

    # Print message
    if (verbose) {
      t1 <- Sys.time()
      msg <- paste0("Extracting ", k, "-mers from ", chr.name)
      dots <- rep(".", 40 - nchar(msg)) |> paste(collapse = "")
      cat(paste0(msg, dots, if(ncpu > 1) "\n"))
    }

    # For removing case k-mer overlap
    if (rm.overlap.region) {

      # Load the coordinate
      coor[chr.name,
           state = ifelse(expand.to.kmer, "kmer", "case"),
           k = ifelse(expand.to.kmer, k, NA)]

      # Mark the coordinate
      coor$mark_overlap()

      # Print message
      if (verbose) {
        percent.overlap <-
          coor[chr.name][, round(sum(is_overlap) / .N * 100, 2)]
        msg2 <- paste0("Removing ", percent.overlap, "% overlaps")
        if (ncpu > 1) msg2 <- paste0(msg2, " in ", chr.name)
        dots2 <- rep(".", 40 - nchar(msg2)) |> paste0(collapse = "")
        cat(paste0(if(ncpu == 1) "\n", msg2, dots2, "\n"))

      }
    }

  # coor[chr.name, state = ifelse(expand.to.kmer, "kmer", "case"), k = ifelse(expand.to.kmer, k, NA)]
  #   map.kmers(chr.name, start, k, end = if(end.coor.matters) end, rm.trunc.kmer = TRUE)

    coor[chr.name,
         state = ifelse(expand.to.kmer, "kmer", "case"),
         k = ifelse(expand.to.kmer, k, NA)][

           # Filter rows
           if (rm.overlap.region) (!is_overlap) else TRUE,

           {



             neg.strand.matters <- coor$is_strand_sensitive && strand == "-"
             end.coor.matters <- !expand.to.kmer & is.null(coor$single_len)

             if (do.standard.extraction | do.growing.kmer.table) {

               # Map all k-mers
               kmers <- map.kmers(chr.name, start, k,
                                  end = if(end.coor.matters) end,
                                  rm.trunc.kmer = TRUE)

               if (neg.strand.matters) {
                 kmers <- reverseComplement(kmers)
               }

               # Count k-mers
               kmers <- data.table(kmer = kmers)[, .N, by = kmer]
               setkey(kmers, kmer)

             } else if (do.split.kmer.extraction) {

               k.1 <- k %/% 2
               k.2 <- k - k.1

               kmers.1 <- map.kmers(
                 chr.name,
                 start = start,
                 k = if(neg.strand.matters) k.2 else k.1,
                 end = if(end.coor.matters) end - k.2,
                 rm.trunc.kmer = FALSE)

               kmers.2 <- map.kmers(
                 chr.name,
                 start = start + if(neg.strand.matters) k.2 else k.1,
                 k = if(neg.strand.matters) k.1 else k.2,
                 end = if(end.coor.matters) end,
                 rm.trunc.kmer = FALSE)

               if (neg.strand.matters) {
                 kmers.1 <- reverseComplement(kmers.1)
                 kmers.2 <- reverseComplement(kmers.2)
               }

               kmers <- data.table(
                 kmer_part1 = if(neg.strand.matters) kmers.2 else kmers.1,
                 kmer_part2 = if(neg.strand.matters) kmers.1 else kmers.2)[
                   , .N, by = .(kmer_part1, kmer_part2)]
               setkey(kmers, kmer_part1, kmer_part2)
             }

             # Update kmer.table for existing k-mers in the table
             if (neg.strand.matters) {
               kmer.table[kmers, neg_strand := neg_strand + N]
             } else {
               kmer.table[kmers, pos_strand := pos_strand + N]
             }

             # Grow kmer.table
             if (do.growing.kmer.table) {

               kmer.col <- names(kmer.table)[grep("kmer", names(kmer.table))]

               # Add additional column to kmers
               if (coor$is_strand_sensitive && strand == "-") {
                 setnames(kmers, "N", "neg_strand")
                 kmers[, pos_strand := 0]
               } else {
                 setnames(kmers, "N", "pos_strand")
                 kmers[, neg_strand := 0]
               }

               # Bind new k-mers found
               kmer.table <- rbind(kmer.table, kmers[!kmer.table])
               setkeyv(kmer.table, kmer.col)

               # Remove different central pattern
               if (!is.null(central.pattern)) {
                 flank.size <- (k - nchar(central.pattern)) / 2
                 rgx <- paste0("^.{", flank.size, "}", central.pattern,
                               collapse = ")|(")
                 rgx <- paste0("(", rgx, ")")
                 kmer.table <- kmer.table[kmer %like% rgx]
                 setkeyv(kmer.table, kmer.col)
               }

               assign("kmer.table", kmer.table, envir = env)
             }

             NULL
           }, by = eval(if(coor$is_strand_sensitive) "strand")]

    # Print message
    if (verbose) {
      t <- Sys.time() - t1
      cat(paste0(if(rm.overlap.region | ncpu > 1) paste0(msg, dots),
                 "DONE! -- ", round(t[[1]], 2), " ", attr(t, "units"), "\n"))
    }

    return(kmer.table)
  }, future.seed = NULL)

  # Aggregate k-mer table if running in different cpu cores
  if (ncpu > 1) {
    for (kmer.count in kmer.counts) {

      if (do.growing.kmer.table) {

        kmer.table <- rbind(kmer.table, kmer.count[!kmer.table])
        setkey(kmer.table, kmer)

      } else {
        kmer.table[kmer.count, `:=`(pos_strand = pos_strand + i.pos_strand,
                                    neg_strand = neg_strand + i.neg_strand)]
      }
    }
  }

  if (!coor$is_strand_sensitive) {

    kmer.table <- countRevCompKmers(kmer.table)

    # If has pattern, remove complementary central.pattern (It was needed for
    # countRevCompKmers operation)
    if (!is.null(central.pattern)) {

      flank.size <- (k - nchar(central.pattern)) / 2
      rgx <- paste0("^.{", flank.size, "}", central.pattern, collapse = ")|(")
      rgx <- paste0("(", rgx, ")")

      kmer.table <- kmer.table[kmer %like% rgx]
    }
  }

  return(kmer.table)
}