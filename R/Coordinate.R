#' Class constructor - build Coordinate object
#'
Coordinate <- R6::R6Class(

  classname = "Coordinate",

  public = list(

    #' @field root_path A path to a directory containing coordinate files.
    root_path = NULL,

    #' @field single_len Single case length e.g. damage length. Default is NULL.
    single_len = NULL,

    #' @field is_strand_sensitive Coordinate strand polarity. Default is TRUE.
    is_strand_sensitive = TRUE,

    #' @field merge_replicate Merge coordinate from different replicates.
    #'     Default is TRUE.
    merge_replicate = TRUE,

    #' @field rm_dup Remove duplicate entry in the coordinate table.
    #'     Default is TRUE.
    rm_dup = TRUE,

    #' @field add_col_rep If add_col_rep is TRUE, column replicate is added to
    #'     the coordinate table. Default is TRUE.
    add_col_rep = FALSE,

    #' @field paths Individual coordinate files.
    paths = NULL,

    #' @field rep_names Replicate names determined from coordinate subdirectory.
    rep_names = NULL,

    #' @field chr_names Chromosome names determined from filenames.
    chr_names = NULL,

    #' @field coor Chromosome-named list of coordinate data.table.
    coor = list(),

    #' @field is_kmer A data.table of is_kmer status. The first column is
    #'     original is_kmer status.
    is_kmer = data.table(original = FALSE),

    #' @field k K-mer size when is_kmer is TRUE. When is_kmer is FALSE, k is NA.
    k = data.table(original = NA),

    #' @field ori_first_index Original chromosome-separated table first index is
    #'     either starting from zero or one.
    ori_first_index = 1,

    #' @field load_limit Maximum coordinate table loaded.
    load_limit = 1,

    #' @description
    #' Create a new Coordinate class
    #' @param root.path A path to a directory containing either:
    #'     (1) chromosome-separated coordinate files
    #'         (assume replicates for subdirectories) OR
    #'     (2) bedfile. (assume replicates for bedfiles)
    #' @param single.len Single case length e.g. damage length. Default is NULL
    #' @param is.kmer Is the coordinate refers to k-mer i.e. expanded case?
    #'     Default is FALSE.
    #' @param k Length of k-mer if is_kmer is TRUE.
    #' @param is.strand.sensitive A boolean whether strand polarity matters.
    #'     Default is TRUE.
    #' @param merge.replicate Merge coordinate from different replicates.
    #'     Default is TRUE. If not merging, duplicates will give weight to the
    #'     kmer counting. If add_col_rep, merged coordinate will contain
    #'     column replicate e.g. "rep1&rep2".
    #' @param add.col.rep Add column replicate to coordinate table.
    #' @param rm.dup Remove duplicates in each replicate. Default is FALSE
    #'     Default is FALSE
    #' @param load.limit Maximum coordinate data.table loaded. Default is 1.
    #' @param ori.first.index Zero- or one-based index. Default is 1.
    #' @return A new `Coordinate` object.
    initialize =
      function(root.path, single.len, is.strand.sensitive, merge.replicate,
               rm.dup, add.col.rep, is.kmer, k, ori.first.index, load.limit) {

        if (!missing(root.path)) self$root_path <- root.path
        if (!missing(single.len)) self$single_len <- single.len
        if (!missing(is.strand.sensitive))
          self$is_strand_sensitive <- is.strand.sensitive
        if (!missing(merge.replicate)) self$merge_replicate <- merge.replicate
        if (!missing(rm.dup)) self$rm_dup <- rm.dup
        if (!missing(add.col.rep)) self$add_col_rep <- add.col.rep
        if (!missing(is.kmer)) self$is_kmer <- data.table(original = is.kmer)
        if (!missing(k))
          self$k <- data.table(original = ifelse(self$is_kmer$original, k, NA))
        if (!missing(ori.first.index)) self$ori_first_index <- ori.first.index
        if (!missing(load.limit)) self$load_limit <- load.limit

        # Error checking
        stopifnot(is.character(self$root_path))

        # If not folder but a direct bedfile path
        bedfile <- grepl(private$bed_pattern, self$root_path)
        if (bedfile) {
          self$paths <- self$root_path
          private$convert_bed()
        }

        # Convert bedfile to chromosome-separated files in a folder.
        bedfile <- all(grepl(private$bed_pattern, list.files(self$root_path)))
        if (bedfile) {
          self$paths <- list.files(self$root_path, self$bed_pattern,
                                   full.names = TRUE)
          private$convert_bed()
        }

        self$paths <- private$get_coor_paths()
        self$rep_names <- unique(sub("^.*/", "", dirname(self$paths)))

        self$chr_names <- private$detect_chr_names()

      },

    #' @description
    #' Calling coordinate table by loading on demand. Maximum load is determine
    #'     by load_limit field.
    #' @param chr.name Chromosome name. It can be a vector of chromosomes.
    #' @param state Coordinate state: "current", "case", "kmer". The coordinate
    #'     state is changed automatically on demand. Default is "current".
    #' @param k K-mer size. If state is "kmer", k is needed to expand the
    #'     coordinate.
    #' @param reload Reload the coordinate table from the root.path.
    #'     Default is TRUE.
    #' @param rm.other.cols Remove unnecessary columns for kmertone operation.
    #' @return A single or list of data.table coordinate of requested
    #'     chromosome.
    `[` = function(chr.name, state="current", k, reload=FALSE,
                   rm.other.cols=TRUE) {

      # Check for number of chr.name asked
      if (length(chr.name) > self$load_limit) {
        stop("Requested chromosome coordinate tables exceed load_limit.")
      }

      if (any(!chr.name %in% self$chr_names)) {

        missing.chrs <- chr.name[!chr.name %in% self$chr_names]
        stop("No chromosome name ", paste(missing.chrs, collapse = " "))

      } else if (reload) {

        private$trim_loads(chr.name)
        self$coor[chr.name] <- private$load_coor(chr.name, rm.other.cols)

      } else if (any(!chr.name %in% names(self$coor))) {

        private$trim_loads(chr.name)
        new.chrs <- chr.name[!chr.name %in% names(self$coor)]

        self$coor[new.chrs] <- private$load_coor(new.chrs, rm.other.cols)

      }

      # Expand or contract to k size if applicable
      if (state == "case") {

        is.kmer <- private$get_kmer_status(chr.name)$is.kmer
        chr.in.kmer <- chr.name[chr.name %in% chr.name[is.kmer]]

        if (length(chr.in.kmer) > 0) {
          exist.k <- private$get_kmer_status(chr.in.kmer)$k
          private$kmerize(chr.in.kmer, exist.k, revert = TRUE)
        }

      } else if (state == "kmer") {

        if (missing(k)) stop("k-mer state is requested without input k.")

        is.kmer <- private$get_kmer_status(chr.name)$is.kmer
        chr.in.kmer <- chr.name[chr.name %in% chr.name[is.kmer]]

        if (length(chr.in.kmer) > 0) {
          private$kmerize(chr.in.kmer, revert = TRUE)
        }

        private$kmerize(chr.name, k)

      }

      if (length(chr.name) > 1) {
        return(self$coor[chr.name])
      } else {
        return(self$coor[[chr.name]])
      }
    },

    #' @description
    #' Mark overlapping regions in the coordinate table. A column name
    #'      is_overlap is added.
    #' @param chr.names Chromosome names
    #' @return New column is_overlap is added.
    mark_overlap = function() {

      for (chr.name in names(self$coor)) {

        setorderv(self$coor[[chr.name]],
                  c("start", if(is.null(self$single_len)) "end"))

        self$coor[[chr.name]][
          ,
          group := cumsum(c(1, cummax(head(

            # end coordinate
            if(self$is_kmer[[chr.name]]) {
            end <- start + self$k[[chr.name]] - 1
          } else if (!is.null(self$single_len)) {
            end <- start + self$single_len - 1
          } else {
            end
          },

          -1)) - tail(start, -1) < 0))
        ]

        self$coor[[chr.name]][, is_overlap := .N > 1, by = group]

        self$coor[[chr.name]][, group := NULL]

      }

    },

    #' @description
    #' Print `Coordinate` object parameters.
    #' @return Message of `Coordinate` object parameters.
    print = function() {

      attrs <- data.table(
        Attribute = c("Root path", "Total replicate", "Case length",
                      "Strand sensitive", "Merge replicate", "Remove duplicate",
                      "Has column replicate", "Total chromosome",
                      "Max table load"),
        Value = c(self$root_path,
                  length(self$rep_names),
                  if(is.null(self$single_len)) "Varied" else self$single_len,
                  self$is_strand_sensitive,
                  self$merge_replicate,
                  self$rm_dup,
                  self$add_col_rep,
                  length(self$chr_names),
                  self$load_limit)
      )

      cat("Coordinate information:\n")
      print(attrs, row.names = FALSE, justify = "left", col.names = "none")

      setcolorder(self$is_kmer, c("original", names(self$coor)))
      cat("\nCoordinate is in k-mer state?\n")
      print(self$is_kmer, row.names = FALSE, justify = "left")

      setcolorder(self$k, c("original", names(self$coor)))
      cat("\nCorresponding k-mer length:\n")
      print(self$k, row.names = FALSE, justify = "left")

      cat("\nLoaded chromosome case:", if(length(self$coor) == 0) "None" else
        names(self$coor), "\n")
      if(length(self$coor) > 0) print(self$coor)

    },

    #' @description
    #' Get corresponding sequence from the loaded coordinate.
    #' @param genome Genome object or vector of named chromosome sequences.
    #' @return New column seq.
    map_sequence = function(genome) {

      for (chr.name in names(self$coor)) {
        self$coor[[chr.name]][,
          seq := stri_sub(genome[chr.name], start,
                          if (is.null(self$single_len)) end
                          else if (self$is_kmer[[chr.name]])
                            start + self$k[[chr.name]] - 1
                          else start + self$single_len - 1)]
        if (self$is_strand_sensitive)
          self$coor[[chr.name]][strand == "-", seq := reverseComplement(seq)]
      }
    }
  ),

  private = list(

    # @field coor_file_pattern. A coordinate file extension pattern.
    coor_file_pattern = "\\.(csv|txt|tsv)($|\\.gz$)",

    # @field bed_pattern. A bed extension pattern.
    bed_pattern = "\\.bed($|\\.gz)",

    # @field chr_pattern A chromosome file pattern.
    chr_pattern = "chr([0-9]{1,2}|[IXV]{1,4}|[XYM])",

    # @description
    # Convert bedfiles to chromosome-separated files.
    # @return None
    convert_bed = function() {

      new.root.path <- tempfile()

      for (rep in seq_along(self$paths)) {
        new.path <- paste0(new.root.path, "/",
                           if(length(self$paths) > 1) paste0("rep_", rep, "/"))
        dir.create(new.path, recursive = TRUE, showWarnings = FALSE)
        bedToCoor(self$paths[rep], new.path)
      }
        self$root_path <- new.root.path
        self$ori_first_index <- 1
    },

    # @description
    # Get coordinate files with full name.
    # @param chr.name Chromosome name. It can be a vector of chromosomes.
    # @return A vector or single coordinate files with full name.
    get_coor_paths = function(chr.name=NULL) {

      if (is.null(chr.name)) {
        chr.names <- private$chr_pattern
      } else {
        chr.names <- paste(chr.name, collapse = "|")
      }

      list.files(self$root_path,
                 paste0("(", chr.names, ")", private$coor_file_pattern),
                 full.names = TRUE, recursive = TRUE)

    },

    # @description
    # Detect chromosome names based on filenames.
    # @return A vector of chromosome names.
    detect_chr_names = function() {
      chr.names <- stri_extract_first_regex(basename(self$paths),
                                            private$chr_pattern) |> unique()

      # Sort numeric first
      chr.no <- sub("chr", "", chr.names) |>
        stri_extract_first_regex("^[0-9]+$") |>
        as.numeric() |> na.omit() |> sort()

      # If zero length, assume in roman
      if (length(chr.no) == 0) {
        chr.no <- sub("chr", "", chr.names) |>
          stri_extract_first_regex("^[IXV]+$") |>
          as.roman() |> na.omit() |> sort()
      }

      # Add back chr prefix
      if (length(chr.no) > 0) {
        chr.no <- paste0("chr", chr.no)
      }

      # Combine all numeric and non-numeric chromosome names
      chr.names <- c(chr.no, chr.names[!chr.names %in% chr.no])

      return(chr.names)
    },

    # @description
    # Update coordinate k-mer status
    # @param chr.name Chromosome name. It can be a vector of chromosomes.
    # @param is_kmer A new boolean for is_kmer parameter.
    # @param k A new k-mer size. k is set to NA when is_kmer is FALSE.
    # @return None
    update_kmer_status = function(chr.name, is.kmer, k) {

      self$is_kmer[, (chr.name) := is.kmer]
      self$k[, (chr.name) := if(!is.kmer) NA else k]

    },

    get_kmer_status = function(chr.name) {

      is.kmer <- unlist(self$is_kmer[, ..chr.name])
      k <- unlist(self$k[, ..chr.name])

      return(list(is.kmer = is.kmer, k = k))
    },

    load_coor = function(chr.name, rm.other.cols=TRUE) {

      coor.paths <- private$get_coor_paths(chr.name)

      coor <- lapply(coor.paths, function(coor.path) {

        # data.table::fread segfault so use nThread=1
        coor <- fread(coor.path, showProgress = FALSE, nThread = 1)

        sample.num <- min(100, nrow(coor))

        # Resolve start coordinate
        if (any(idx <- grepl("^start|^Start|^START", names(coor)))) {
          setnames(coor, names(coor)[idx], "start")
        } else {
         # expect the first column of round numbers is start
          col.start <- coor[1:sample.num][,
            names(.SD)[sapply(.SD, function(e) {
              if (is.numeric(e)) all(e %% 1 == 0) else FALSE
            })]][1]
          setnames(coor, col.start, "start")
        }

        # Resolve chromosome name
        chr.name <- stri_extract_first_regex(basename(coor.path),
                                             "chr[0-9A-Za-z]+")

        # Resolve replicate name if any (sub-folder)
        rep.name <- sub("^.*/", "", dirname(coor.path))
        if (rep.name == coor.path) rep.name <- NULL

        # Resolve end coordinate
        if (is.null(self$single_len)) { # expect column end

          if (any(idx <- grepl("^end|^End|^END", names(coor)))) {

            setnames(coor, names(coor)[idx], "end")

          } else {
            # expect the second column of round numbers is end
             col.end <- coor[1:sample.num][, names(.SD)[
                 sapply(.SD, function(e) {
                   if (is.numeric(e)) all(e %% 1 == 0) else FALSE
                 })]][2]
            if (is.na(col.end)) {
              stop("Failed to find column end. Please check the coordinate ",
                   "table and single.len argument.")
            }
            setnames(coor, col.end, "end")
          }

        } else {
          if ("end" %in% names(coor)) coor[, end := NULL]
        }

        # Resolve strand polarity
        if (self$is_strand_sensitive) {
          if (any(idx <- grepl("^strand|^Strand|^STRAND", names(coor)))) {
            setnames(coor, names(coor)[idx], "strand")
          } else if (any(idx <- coor[1:sample.num][
            , sapply(.SD, function(col.elm) {
              all(col.elm %in% c("+", "-", "*"))
            })])) {
            setnames(coor, names(coor)[idx], "strand")
          } else {
            stop("Coordinate is.strand.sensitive but there is no strand column")
          }
        }

        # Resolve necessary columns
        col.names <- c("start",
                       if (is.null(self$single_len)) "end",
                       if (self$is_strand_sensitive) "strand")

        if (rm.other.cols) {
          for (col.name in names(coor)[!names(coor) %in% col.names]) {
            set(coor, j = col.name, value = NULL)
          }
        }
        setcolorder(coor, col.names)

        # Remove duplicates
        if (self$rm_dup) {
          coor <- unique(coor)
        }
        setorderv(coor, col.names)
        setkeyv(coor, col.names)

        # Return useful contents
        coor <- list(coor = coor, chr.name = chr.name, rep.name = rep.name)

        return(coor)
      })

      # Separate rep.names and coor
      chr.names <- sapply(coor, `[[`, "chr.name")
      rep.names <- unique(unlist(lapply(coor, `[[`, "rep.name")))
      coor <- lapply(coor, `[[`, "coor")

      names(coor) <- chr.names
      chr.names <- unique(chr.names)

      # Merge or combine replicates if any
      if (!is.null(rep.names)) {

        coor <- lapply(chr.names, function(chr.name) {

          if (self$add_col_rep) {
            for (i in seq_along(rep.names)) {
              coor[names(coor) == chr.name][[i]][
                , paste0("rep", i) := rep.names[i]]
            }
          }

          if (self$merge_replicate) {

            coor <- Reduce(function(x, y) merge.data.table(x, y, all = TRUE),
                           coor[names(coor) == chr.name])

            if (self$add_col_rep) {

              rep.cols <- names(coor)[grep("rep[0-9]+", names(coor))]

              rep.sep <- "&"
              coor[, replicate := do.call(paste, c(.SD, sep = rep.sep)) |>
                     stri_replace_all_regex(paste0("(NA", rep.sep, ")|(",
                                                   rep.sep, "NA)"),
                                            ""),
                   .SDcols = rep.cols]

              set(coor, j = rep.cols, value = NULL)
            }

          } else {
            coor <- rbindlist(coor[names(coor) == chr.name])
          }

          if (self$ori_first_index == 0) coor[, start := start + 1]

          return(coor)
        })

        names(coor) <- chr.names
      }

      # Set original k-mer status
      self$is_kmer[, (chr.name) := original]
      self$k[, (chr.name) := original]

      return(coor)
    },

    kmerize = function(chr.names, k, revert=FALSE) {

      check.k <- function(flank.size, k) {
        if (flank.size * 2 + self$single_len != k) {
          stop("Unbalanced flank size to form k-mer with the case at the ",
               "center. Please check size k.")
        }
      }

      for (chr.name in chr.names) {

        # Contract to case coordinate
        if (self$is_kmer[[chr.name]] & revert) {

          flank <- (self$k[[chr.name]] - self$single_len) %/% 2
          check.k(flank, self$k[[chr.name]])
          self$coor[[chr.name]][, start := start + flank]
          private$update_kmer_status(chr.name, is.kmer = FALSE, k = NA)

          # Expand to k-mer coordinate
        } else if (!self$is_kmer[[chr.name]] & self$single_len < k & !revert) {

          flank <- (k - self$single_len) %/% 2
          check.k(flank, k)
          self$coor[[chr.name]][, start := start - flank]
          private$update_kmer_status(chr.name, is.kmer = TRUE, k = k)

        } else if (self$single_len > k) {

          stop("The case length is larger than k. Cannot expand to k-mer.")

        }
      }
    },

    trim_loads = function(chr.name) {

      exist.chrs <- names(self$coor)
      new.chrs <- chr.name[!chr.name %in% exist.chrs]

      if (length(self$coor) + length(new.chrs) > self$load_limit) {
        avail.slots <- which(!exist.chrs %in% chr.name)
        self$coor[idx.rm <- avail.slots[length(new.chrs)]] <- NULL
        set(self$is_kmer, j = exist.chrs[idx.rm], value = NULL)
        set(self$k, j = exist.chrs[idx.rm], value = NULL)
      }
    }

  ))

#' @export
`[.Coordinate` <- function(obj, ...) obj$`[`(...)
