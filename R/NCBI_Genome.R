#' Class constructor - build NCBI Genome object
#'
#' NCBI FASTA file contain nucleotide accession number at the headers, followed
#'    by some information about the sequence whether they are chromosome,
#'    plasmid, or mictochondria, their assembly status, etc.
NCBI_Genome <- R6::R6Class(

  classname = "NCBI_Genome",

  public = list(

    #' @field fasta_file A path to FASTA file.
    #'     fasta files.
    fasta_file = NULL,

    #' @field genome_name A genome name.
    genome_name = "unknown",

    #' @field db NCBI database: "refseq" or "genbank"
    db = "unknown",

    #' @field seq A chromosome-named list of sequences.
    seq = list(),

    #' @field seq_len A chromosome-named vector of sequence length.
    seq_len = NULL,

    #' @field load_limit Maximum chromosome sequences loaded.
    load_limit = 1,

    #' @field mask Genome mask status: "hard", "soft", or "none".
    mask = "none",

    #' @field use_UCSC_name Use UCSC style chromosome name? Default to FALSE.
    use_UCSC_name = FALSE,

    #' @field headers A chromosome-named vector of headers.
    headers = NULL,

    #' @field avail_seqs Available chromosome sequences in the fasta file.
    avail_seqs = NULL,

    #' @field asm Assembly summary.
    asm = NULL,

    #' @description
    #' Create a new NCBI Genome class
    #' @param genome.name A genome name. NCBI genome is included with kmeRtone.
    #' @param db NCBI database: "refseq" or "genbank".
    #' @param fasta.file A path to the NCBI-style fasta files. This is for
    #'    user's own FASTA file.
    #' @param asm NCBI assembly summary.
    #' @param load.limit Maximum chromosome sequences loaded. Default is 1.
    #' @param use.UCSC.name Use UCSC style chromosome name? Default to FALSE.
    #' @param mask Genome mask status: "hard", "soft", or "none". Default is
    #'     "none".
    #' @return A new `NCBI Genome` object.
    initialize = function(genome.name, db, fasta.file, asm, mask, use.UCSC.name,
                          load.limit) {

      if (!missing(genome.name)) self$genome_name <- genome.name
      else if (missing(genome.name)) self$genome_name <- basename(fasta.file)
      if (!missing(use.UCSC.name)) self$use_UCSC_name <- use.UCSC.name
      if (!missing(fasta.file)) self$fasta_file <- fasta.file
      if (!missing(load.limit)) self$load_limit <- load.limit
      if (!missing(mask)) self$mask <- mask
      if (!missing(asm)) self$asm <- asm

      # Guessing db if missing
      if (!missing(db)) self$db <- db
      else if (missing(db)) {
        prefix <- stri_sub(genome.name, 1, 4)
        if (prefix == "GCF_") {
          self$db <- "refseq"
        } else if (prefix == "GCA_") {
          self$db <- "genbank"
        } else if (missing(fasta.file)) {
          stop('Please input NCBI database name: "refseq" or "genbank"')
        }
      }

      if (is.null(self$fasta_file)) {
        self$fasta_file <- private$get_fasta_path(self$genome_name, self$db)
      }

      # Resolve chromosome name by building data.table of headers.
      private$dt_headers <- private$build_dt_headers(self$fasta_file,
                                                     self$use_UCSC_name)
      self$headers <- private$dt_headers[-.N, header]
      names(self$headers) <- private$dt_headers[-.N, chromosome]
      self$avail_seqs <- names(self$headers)

    },

    #' @description
    #' Calling chromosome sequence by loading on demand.
    #'     Maximum load is determine by load_limit field.
    #' @param chr.names Chromosome name. It can be a vector of chromosomes.
    #' @param reload Reload the sequence from the fasta_file.
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
        self$seq_len <- self$seq_len[names(self$seq_len) %in% names(self$seq)]

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

      cat("NCBI Genome:", self$genome_name, "\n")
      cat("Available sequence:", if(length(self$avail_seqs) < 30)
        self$avail_seqs else c(self$avail_seqs[1:30], "....."), "\n")
      cat("Loaded sequence:", if(length(self$seq) == 0) "None\n" else "\n")
      if(length(self$seq) > 0) print(self$seq_len)
    },

    #' @description
    #' Get NCBI assembly report for the genome.
    #' @return Message of `Genome` object summary.
    get_assembly_report = function() {
      report.path <- gsub("_genomic.fna.gz", "_assembly_report.txt",
                          self$fasta_file)
      fread(report.path, skip = "Sequence-Name", showProgress = FALSE)
    }),

  private = list(

    # @field dt_headers FASTA file headers with their corresponding line number.
    dt_headers = data.table(),

    # @description
    # Build metadata data.table of header for skipping row when reading fasta.
    # @return A data.table of headers with chromosome name and line number in
    #    the fasta file.
    build_dt_headers = function(fasta.file, use.UCSC.name) {

      index.file <- gsub(".fna.gz", ".index.gz", self$fasta_file)

      # Check the index is already built.
      if (file.exists(index.file)) {
        dt.headers <- fread(index.file, showProgress = FALSE)
        return(dt.headers)
      }

      # Growing row
      dt.headers <- data.table(chromosome = character(), header = character(),
                               line_number = numeric())
      # Growing column
      # dt.headers <- data.table(column.name = c("chromosome", "header",
      #                                          "line_number"))
      #cat("Indexing the fasta file...\n")
      # Assume the genome is large, so process by chunk.
      max.row <- 12000000
      skip.len <- 0
      chunk.num <- 0
      fasta.file.conn <- file(fasta.file, "rt")
      while (TRUE) {
        fasta <- scan(fasta.file.conn, what = "character", sep = "\n",
                      nlines = max.row, quiet = TRUE, quote = "")
        line.numbers <- stri_startswith_fixed(fasta, ">") |> which()
        if (length(line.numbers) == 0) {
          chunk.num <- chunk.num + 1
          skip.len <- chunk.num * max.row
          next
        }
        headers <- fasta[line.numbers]
        line.numbers <- line.numbers + skip.len
        chr.names <- gsub(">", "", stri_extract_first_regex(headers, ">\\S+"))
        # Growing column
        # for (i in seq_along(chr.names))
        #   set(x = dt.headers, j = paste0(chunk.num, "_", i),
        #       value = c(chr.names[i], headers[i], line.numbers[i]))
        # Growing row
        dt.headers <- rbind(dt.headers,
                            data.table(chromosome = chr.names,
                                       header = headers,
                                       line_number = line.numbers))

        if (length(fasta) < max.row) break
        chunk.num <- chunk.num + 1
        skip.len <- chunk.num * max.row
      }
      close(fasta.file.conn)
      # Growing column
      # dt.headers[, endline := c(NA, NA, skip.len + length(fasta) + 1)]
      # dt.headers <- transpose(dt.headers, make.names = "column.name")
      # dt.headers[, line_number := as.numeric(line_number)]
      # Growing row
      dt.headers <- rbind(dt.headers,
        data.table(chromosome = NA, header = NA, line_number = skip.len +
          length(fasta) + 1))
      if (use.UCSC.name) {
        report <- self$get_assembly_report()
        setnames(report, names(report), stri_trans_tolower(names(report)))
        setkeyv(report, paste0(db, "-accn"))
        setkey(dt.headers, chromosome)
        dt.headers[report, chromosome := `ucsc-style-name`]
      }
      setkey(dt.headers, line_number)
      fwrite(dt.headers, gsub(".fna.gz", ".index.gz", self$fasta_file),
             showProgress = FALSE)
      return(dt.headers)
    },

    # @description
    # Get FASTA path for a given genome name. If not found, try to download
    #    from NCBI server.
    # @param genome.name Genome name.
    # @param db NCBI database: "refseq" or "genbank"
    # @return fasta.file path
    get_fasta_path = function(genome.name, db) {
      genome.dir <- paste0(path.expand(kmeRtone.data.path), "/genome/", db)
      fasta.file <- list.files(genome.dir, genome.name,
                               full.names = TRUE)
      fasta.file <- fasta.file[grepl("_genomic.fna.gz", fasta.file)]
      if (length(fasta.file) > 1) {
        cat("More than one fasta files are found:\n", basename(fasta.file),
            "\nThis version is used:", basename(fasta.file[1]), "\n")
        fasta.file <- fasta.file[1]
      } else if (length(fasta.file) == 0) {
        cat(genome.name, "is not found in local database. Looking up in remote",
            db, "database...\n")
        if (is.null(self$asm)) {
          self$asm <- getNCBIassemblySummary(organism.group = "all", db = db)
        }
        self$asm <- self$asm[ftp_path %like% genome.name] |>
          selectRepresentativeFromASM()
        if (nrow(self$asm) == 0) {
          stop(genome.name, " not found in ", db, " database.")
        } else {
          cat(genome.name, "found! Downloading...\n")
          downloadNCBIGenomes(self$asm, output.dir = genome.dir)
        }
        fasta.file <- private$get_fasta_path(self$genome_name, self$db)
      }
      return(fasta.file)
    },

    # @description
    # Load chromosome sequence on demand.
    # @param chr.name A single or vector of chromosome names.
    # @param mask mask i.e. capitalize the sequence. Only applicable to
    #    soft-masked genome. UCSC genomes are soft masked.
    # @return A list of sequence(s).
    load_seq = function(chr.names, mask=self$mask) {
      seq <- lapply(chr.names, function(chr.name) {
        idx <- which(private$dt_headers$chromosome == chr.name)
        if (length(idx) == 0) stop(chr.name, " not found!")
        seq <- scan(self$fasta_file, what = "character", sep = "\n",
                    skip = private$dt_headers[idx, line_number],
                    nlines = private$dt_headers[c(idx, idx + 1),
                                                abs(diff(line_number)) - 1],
                    quiet = TRUE, quote = "") |> stri_paste(collapse = "")
        if (mask == "none") {
          seq <- stri_trans_toupper(seq)
        } else if (mask == "hard") {
          seq <- stri_sub_all_replace(seq, stri_locate_all_regex(seq, "[a-z]"),
                                      replacement = "N")
        }
      })
      return(seq)
    },

    # @description
    # A utility function to remove chromosome when maximum load reaches.
    # @param chr.name A single or vector of newly requested chromosome names.
    # @return None.
    trim_loads = function(chr.name) {

      new.chrs <- chr.name[!chr.name %in% names(self$seq)]

      if (length(self$seq) + length(new.chrs) > self$load_limit) {
        avail.slots <- which(!names(self$seq) %in% chr.name)
        self$seq[avail.slots[length(new.chrs)]] <- NULL
      }
    }

  ))

#' @export
`[.NCBI_Genome` <- function(obj, ...) obj$`[`(...)
