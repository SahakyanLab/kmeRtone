#' Class constructor - build Genome object
#' @importFrom R6 R6Class
#' @importFrom data.table fread fwrite setorder 
#' @importFrom stringi stri_length stri_count_fixed stri_sub
UCSC_Genome <- R6::R6Class(

  classname = "UCSC_Genome",

  public = list(

    #' @field root_path A path to a directory containing chromosome-separated
    #'     fasta files.
    root_path = NULL,

    #' @field genome_name A genome name.
    genome_name = "unknown",

    #' @field paths Individual chromosome sequence files.
    paths = NULL,

    #' @field seq A chromosome-named list of sequences.
    seq = list(),

    #' @field seq_len A chromosome-named vector of sequence length.
    seq_len = NULL,

    #' @field load_limit Maximum chromosome sequences loaded.
    load_limit = 1,

    #' @field mask Genome mask status: "hard", "soft", or "none".
    mask = "none",

    #' @field info_file Path to info file with pre-computed values.
    info_file = "info.csv",

    #' @field chr_names Chromosome names.
    chr_names = NULL,

    #' @description
    #' Create a new Genome class
    #' @param genome.name A genome name. UCSC genome is included with kmeRtone.
    #' @param root.path A path to a directory containing chromosome-separated
    #'     fasta files.
    #' @param load.limit Maximum chromosome sequences loaded. Default is 1.
    #' @param mask Genome mask status: "hard", "soft", or "none". Default is
    #'     "none".
    #' @return A new `Genome` object.
    initialize = function(genome.name, root.path, mask, load.limit) {

      if (!missing(genome.name)) self$genome_name <- genome.name
      else if (missing(genome.name)) self$genome_name <- basename(root.path)
      if (!missing(root.path)) self$root_path <- root.path
      if (!missing(load.limit)) self$load_limit <- load.limit
      if (!missing(mask)) self$mask <- mask

      if (is.null(self$root_path)) {
        home.path <- path.expand(kmeRtone.data.path)
        self$root_path <- paste0(home.path, "/genome/", self$genome_name, "/")
      }

      self$paths <- private$get_seq_path()
      private$avail_chrs <- private$detect_local_chrs()

      self$info_file <- paste0(self$root_path, "/", self$info_file)

      if (file.exists(self$info_file)) {
        self$chr_names <- fread(self$info_file, showProgress = FALSE)$chromosome
      }

      # Error checking
      if (length(self$paths) == 0 &
          !grepl(private$available_genomes, self$genome_name)) {
        stop("No fasta files detected. Please check your genome root path",
             " or genome name.")
      }
    },

    #' @description
    #' Calling chromosome sequence by loading on demand.
    #'     Maximum load is determine by load_limit field.
    #' @param chr.names Chromosome name. It can be a vector of chromosomes.
    #' @param reload Reload the sequence from the root_path.
    #'     Default is FALSE.
    #' @return A single or list of sequence of requested chromosome.
    `[` = function(chr.names, reload=FALSE) {

      # Check for number of chr.names asked
      if (length(chr.names) > self$load_limit) {
        stop("Requested chromosome data exceed load_limit.")
      }

      private$trim_loads(chr.names)

      if (reload) {

        self$seq[chr.names] <- private$load_seq(chr.names, self$mask)
        self$seq_len[chr.names] <- stri_length(self$seq[chr.names])

      } else if (any(!chr.names %in% names(self$seq))) {

        new.chrs <- chr.names[!chr.names %in% names(self$seq)]
        self$seq[new.chrs] <- private$load_seq(new.chrs, self$mask)
        self$seq_len[new.chrs] <- stri_length(self$seq[new.chrs])

      }

     # if (length(chr.name) > 1) {
     #   return(self$seq[chr.name])
     # } else {
     #   return(self$seq[[chr.name]])
     # }
      return(self$seq[chr.names])
    },

    #' @description
    #' Print summary of `Genome` object.
    #' @return Message of `Genome` object summary.
    print = function() {

      chr.nums <- stri_sub(names(self$seq), 4)
      idx <- grepl("[0-9]", chr.nums)
      loaded.chrs <- paste0("chr", c(sort(as.numeric(chr.nums[idx])),
                                     sort(chr.nums[!idx])))

      cat("Genome:", self$genome_name, "\n")
      cat("Available sequence:", private$avail_chrs, "\n")
      cat("Loaded sequence:", if(length(chr.nums) == 0) "None\n" else "\n")
      if(length(self$seq) > 0) {
        print(self$seq_len[loaded.chrs])
      }
    },

    #' @description
    #' Get chromosome length from pre-calculated length
    #' @param chr.names Chromosome name. It can be a vector of chromosomes.
    #' @param recalculate Recalculate the pre-calculated length.
    #' @return A chromosome-named vector of length value.
    get_length = function(chr.names, recalculate=FALSE) {

      info.dt <- private$read_info()

      if (!"length" %in% names(info.dt) | recalculate) {
        info.dt[chromosome %in% chr.names,
                length := stri_length(private$load_seq(chromosome)),
                by = chromosome]
        fwrite(info.dt, self$info_file)
      } else if (info.dt[chromosome %in% chr.names & is.na(length), .N] > 0) {
        info.dt[chromosome %in% chr.names & is.na(length),
                length := stri_length(private$load_seq(chromosome)),
                by = chromosome]
        fwrite(info.dt, self$info_file)
      }

      len <- info.dt$length
      names(len) <- info.dt$chromosome

      return(len[chr.names])
    },

    #' @description
    #' Get pre-calculated sequence content e.g. G+C content
    #' @param chr.names Chromosome name. It can be a vector of chromosomes.
    #' @param seq Sequence to count. e.g. c("G", "C")
    #' @param recalculate Recalculate the pre-calculated length.
    #' @return A chromosome-named vector of sequence content.
    get_content = function(chr.names, seq, recalculate=FALSE) {

      info.dt <- private$read_info()

      if (all(!seq %in% names(info.dt)) | recalculate) {
        info.dt[
          chromosome %in% chr.names,
          eval(seq) := stri_count_fixed(private$load_seq(chromosome), seq,
                                        overlap = TRUE) |> as.list(),
          by = chromosome
        ]
        fwrite(info.dt, self$info_file)

      } else {

        # Create missing column
        if (any(!seq %in% names(info.dt))) {
          info.dt[, eval(seq[!seq %in% names(info.dt)]) := NA]
        }

        # Calculate missing value
        info.dt[chromosome %in% chr.names &
                sapply(seq, function(s) is.na(get(s))) |> rowSums() > 0,

                eval(seq) := stri_count_fixed(private$load_seq(chromosome), seq,
                                              overlap = TRUE) |> as.list(),
              by = chromosome]

        fwrite(info.dt, self$info_file)
      }

      cnt <- info.dt[, rowSums(.SD), .SDcols = seq]
      names(cnt) <- info.dt$chromosome

      return(cnt[chr.names])
    }

  ),

  private = list(

    # @field fasta_pattern. A csv extension pattern.
    fasta_pattern = "\\.(fa|fna|fasta)($|\\.gz$)",

    # @field available_genomes Available genomes in kmeRtone package)
    #     i.e. UCSC database.
    available_genomes = "^(hg|mm|sacCer)[0-9IVX]+",

    # @field avail_chrs Available chromosome sequences in root_path. UCSC genome
    #    sequences are downloaded once on demand.
    avail_chrs = NULL,

    # @description
    # A utility function to detect locally avail_chrs.
    # @return A vector of locally available chromosomes.
    detect_local_chrs = function() {
      sub(private$fasta_pattern, "",
          basename(private$get_seq_path()))
    },

    # @description
    # A utility function to read info file in the genome root directory.
    # @return A data.table of info.
    read_info = function() {

      if (file.exists(self$info_file)) info.dt <- fread(self$info_file)
      else info.dt <- data.table(chromosome = private$avail_chrs)
      return(info.dt)
    },

    # @description
    # A utility function to download UCSC genome.
    # @param chr.name Chromosome name.
    # @return None.
    download_genome = function(chr.name) {

      downloadUCSCgenome(self$genome_name, self$root_path, chr.name,
                         method = "curl")

      if (is.null(self$chr_names)) {
        self$chr_names <- fread(self$info_file, showProgress = FALSE)$chromosome
      }

      # update available chromosome
      private$avail_chrs <- private$detect_local_chrs()
    },

    # @description
    # Load chromosome sequence on demand.
    # @param chr.names A single or vector of chromosome names.
    # @param mask Genome mask status: "hard", "soft", "none"
    # @return A list of sequence(s).
    load_seq = function(chr.names, mask = self$mask) {

      seq.paths <- private$get_seq_path(chr.names)

      if (length(seq.paths) < length(chr.names)) {

        new.chrs <- chr.names[!chr.names %in% private$avail_chrs]

        if (grepl(private$available_genomes, self$genome_name)) {

          private$download_genome(new.chrs)
          self$paths <- private$get_seq_path()
          seq.paths <- private$get_seq_path(chr.names)

        } else {

          stop("No ", new.chrs, " detected. Please check your fasta file name.")

        }
      }

      seq <- lapply(seq.paths, readSingleFASTA, mask = mask)

      return(seq)
    },

    # @description
    # Get individual chromosome sequence path.
    # @param chr.names A single or vector of chromosome names.
    # @return A single or vector of chromosome fasta path.
    get_seq_path = function(chr.names=NULL) {

      if (is.null(chr.names)) {
        chr.names <- "chr.+"
      } else {
        chr.names <- paste(chr.names, collapse = "|")
      }

      list.files(self$root_path,
                 paste0("(", chr.names, ")", private$fasta_pattern),
                 full.names = TRUE)

    },

    # @description
    # A utility function to remove chromosome when maximum load reaches.
    # @param chr.name A single or vector of newly requested chromosome names.
    # @return None.
    trim_loads = function(chr.names) {

      new.chrs <- chr.names[!chr.names %in% names(self$seq)]

      if (length(self$seq) + length(new.chrs) > self$load_limit) {
        avail.slots <- which(!names(self$seq) %in% chr.names)
        self$seq[avail.slots[length(new.chrs)]] <- NULL
      }
    }

  ))

#' @export
`[.UCSC_Genome` <- function(obj, ...) obj$`[`(...)
