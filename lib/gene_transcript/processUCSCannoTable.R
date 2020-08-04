processUCSCannoTable <- function(genome.annotation, type = "refseq",
                                 element = "all", upstream = NULL,
                                 downstream = NULL, genome = NULL) {
  # Process raw genome annotation table from UCSC Genome. Output genomic
  # coordinate of requested element(s)
  # If genome is not supplied, trimming will not be performed.
  #
  # Flags
  # element    <string>  "all", "intron", "exon", "cds", "utr"
  # upstream    <int>    Length of upstream sequence. It can overlap other gene.
  # downstream  <int>    Length of downstream sequence. It can overlap other
  #                      gene.
  # genome     <genome>  For trimming down upstream and downstream if out of
  #                      range.
  #
  # IDEA version number 2: Make a small function that has two flags. One is the
  # main ranges and 2nd is the unwanted ranges. The function will substract
  # unwanted regions from the main ranges. Maybe it will be cleaner this way.
  #
  # The RefSeq annotation table only explicitly gives coordinate for exons.
  # The other elements need to be inferred based on the information given in
  # the annotation table.
  #
  # All the operations in here are vectorized. If the table is big, expect a
  # spike in memory. Using ncbiRefSeq table and genome hg38, the memory is
  # stable at 4-5 GB. I can utilise data.table package to process by chunk if
  # needed.

  # Original table is zero-based open-end index. The indexing system is changed
  # to follow R's system. Critical information based on UCSC Genome website:
  #    Column        Explaination
  #    bin           Indexing field to speed chromosome range queries. (Only
  #                  relavent to UCSC programme)
  #    name 	     Name of gene (usually transcript_id from GTF)
  #    chrom 	     Reference sequence chromosome or scaffold
  #    strand 	     + or - for strand
  #    txStart 	     Transcription start position (or end position for minus
  #                  strand item)
  #    txEnd 	     Transcription end position (or start position for minus
  #                  strand item)
  #    cdsStart      Coding region start (or end position for minus strand item)
  #    cdsEnd 	     Coding region end (or start position for minus strand item)
  #    exonCount     Number of exons
  #    exonEnds      Exon end positions (or start positions for minus strand
  #                  item)
  #    exonStart     Exon start positions (or end positions for minus strand
  #                  item)
  #    name2 	     Alternate name (e.g. gene_id from GTF)
  #    cdsStartStat  Status of CDS start annotation (none, unknown, incomplete,
  #                  or complete) = ('none','unk','incmpl','cmpl')
  #    cdsEndStat    Status of CDS end annotation (none, unknown, incomplete,
  #                  or complete)
  #    exonFrames    Exon frame {0,1,2}, or -1 if no frame for exon (Related to
  #                  codon. Number represents extra bases (modulus of 3) from
  #                  previous exon block brought to a current exon block.)
  #
  # If cdsStart == cdsEnd, that means non-coding sequence.
  #    - maybe cdsStartStat and cdsEndStat == "none" mean the same thing.
  #      maybe exonFrames == "-1," means the same thing.
  #
  # The final output table follows kmertone genomic coordinate format i.e.
  # columns are chromosome, start, end, and strand.

  # Dependency: data.table, stringi

  if ("all" %in% element) element <- c("exon", "intron", "cds", "utr")

  if (class(genome.annotation)[1] == "character"){
    genome.annotation <- fread(genome.annotation)
  } else if (class(genome.annotation)[1] == "data.frame") {
    setDT(genome.annotation)
  } else if (class(genome.annotation)[1] != "data.table"){
    stop(paste("Please input genome.annotation either as a path to the
               annotation table or <data.table>"))
  }

  if (type == "refseq") {
    columns <- c("bin", "name", "chrom", "strand", "txStart", "txEnd",
                 "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                 "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

    # If no column names yet.
    if (sum(names(genome.annotation) %in% columns) < length(columns)) {
      setnames(genome.annotation, names(genome.annotation), columns)
    }
  }

  if (type == "gencode") {
    columns <- c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart",
                 "cdsEnd", "exonCount", "exonStarts", "exonEnds", "proteinID",
                 "alignID")

    # If no column names yet.
    if (sum(names(genome.annotation) %in% columns) < length(columns)) {
      setnames(genome.annotation, names(genome.annotation), columns)
    }
  }

  # Convert to 1-based indexing (+1 to every start coordinate)
  genome.annotation[, c("txStart", "cdsStart") := .(txStart + 1, cdsStart + 1)]

  # BEGIN
  # Select exons, introns, cds, 5'-UTR, 3'-UTR, upstream, downstream
  genomic.coordinate <- genome.annotation[, {

    # --------------------------------------------------------------------------
    # EXONS

    # Exons coordinates (+1 to convert to 1-based indexing)
    exon.starts <- as.numeric(unlist(stri_split_fixed(exonStarts, ",",
                                                      omit_empty = TRUE))) + 1
    exon.ends <- as.numeric(unlist(stri_split_fixed(exonEnds, ",",
                                                    omit_empty = TRUE)))


    # --------------------------------------------------------------------------
    # For vectorisation purpose

    # Repeat strand, cdsStarts, cdsEnds, txStarts, and txEnds for every exons
    # This is to be used in vectorisation of every genes in the table.
    cds.start <- rep(cdsStart, exonCount)
    cds.end <- rep(cdsEnd, exonCount)
    exon.strand <- rep(strand, exonCount)

    # Exon strand index for + and -
    idx.exon.pos <- exon.strand == "+"
    idx.exon.neg <- exon.strand == "-"


    # --------------------------------------------------------------------------
    # INTRONS

    if("intron" %in% element){

      # First and last exon indices (All exons from different transcripts are
      # combined in a vector)
      idx.last.exons <- cumsum(exonCount)
      idx.first.exons <- shift(idx.last.exons, fill = 0) + 1

      # Introns coordinates
      intron.starts <- exon.ends[-idx.last.exons] + 1
      intron.ends <- exon.starts[-idx.first.exons] - 1

    } else {
      intron.starts <- intron.ends <- NULL
    }

    # --------------------------------------------------------------------------
    # CDS - exons without UTR

    if("cds" %in% element){

      # Indexes of exons that contain cds can be divided into 2 categories:
      # completely cds and partly cds
      ## Upstream indexes of UTR
      idx.contain.cds <- exon.ends >= cds.start & exon.starts <= cds.end
      idx.completely.cds <- exon.starts >= cds.start & exon.ends <= cds.end
                         # complete&partial    complete
                         #       TRUE      +     FALSE
      # idx.partly.cds <- (idx.contain.cds + idx.completely.cds) == 1
      idx.partly.cds.up <- (exon.starts < cds.start & exon.ends > cds.start) |
                           (exon.ends == cds.start)
      idx.partly.cds.down <- (exon.starts < cds.end & exon.ends > cds.end) |
                             (exon.starts == cds.end)

      # In a situation where idx.partly.up and idx.partly.down is the same, cds
      # is in only one exon block
      idx.partly.one <- idx.partly.cds.up & idx.partly.cds.down

      # CDS
      cds.starts <- c(cds.start[idx.partly.cds.up & !idx.partly.one],
                      exon.starts[idx.completely.cds],
                      exon.starts[idx.partly.cds.down & !idx.partly.one],
                      cds.start[idx.partly.one])
      cds.ends   <- c(exon.ends[idx.partly.cds.up & !idx.partly.one],
                      exon.ends[idx.completely.cds],
                      cds.end[idx.partly.cds.down & !idx.partly.one],
                      cds.end[idx.partly.one])

      # Note: cds.start vs. cds.starts (cds.end vs. cds.ends)
      # cds.start is the very start coordinate of CDS. It is the repetition of
      # cdsStart to aid with the vectorisation operation.
      # cds.starts are multiple start coordinates correspond to exon blocks

    } else {
      cds.starts <- cds.ends <- NULL
    }

    # --------------------------------------------------------------------------
    # 5'-UTR, 3'-UTR

    if ("utr" %in% element) {

      # Indexes of exons that contain UTR can be divided into 2 categories:
      # completely UTR and partly UTR
      ## Upstream indexes of UTR
      idx.contain.UTR.up <- exon.starts < cds.start
      idx.completely.UTR.up <- exon.ends < cds.start
                                #      complete&partial complete
                                #    TRUE      +      FALSE
      idx.partly.UTR.up <- (idx.contain.UTR.up + idx.completely.UTR.up) == 1

      ## Downstream indexes of UTR
      idx.contain.UTR.down <- exon.ends > cds.end
      idx.completely.UTR.down <- exon.starts > cds.end
      idx.partly.UTR.down <- (idx.contain.UTR.down +
                              idx.completely.UTR.down) == 1

      # Note: The directionality of 5'-UTR and 3'-UTR are reversed in - strand.
      #       So, the strand information is required here.

      # 5'-UTR coordinates
      UTR5.starts <- c(exon.starts[idx.completely.UTR.up & idx.exon.pos],
                       exon.starts[idx.partly.UTR.up & idx.exon.pos],
                       exon.starts[idx.completely.UTR.down & idx.exon.neg],
                       cds.end[idx.partly.UTR.down & idx.exon.neg] + 1)
      UTR5.ends   <- c(exon.ends[idx.completely.UTR.up & idx.exon.pos],
                       cds.start[idx.partly.UTR.up & idx.exon.pos] - 1,
                       exon.ends[idx.completely.UTR.down & idx.exon.neg],
                       exon.ends[idx.partly.UTR.down & idx.exon.neg])

      # 3'-UTR coordinates
      UTR3.starts <- c(cds.end[idx.partly.UTR.down & idx.exon.pos]+1,
                       exon.starts[idx.completely.UTR.down & idx.exon.pos],
                       exon.starts[idx.partly.UTR.up & idx.exon.neg],
                       exon.starts[idx.completely.UTR.up & idx.exon.neg])
      UTR3.ends   <- c(exon.ends[idx.partly.UTR.down & idx.exon.pos],
                       exon.ends[idx.completely.UTR.down & idx.exon.pos],
                       cds.start[idx.partly.UTR.up & idx.exon.neg] - 1,
                       exon.ends[idx.completely.UTR.up & idx.exon.neg])

    } else {
      UTR5.starts <- UTR5.ends <- UTR3.starts <- UTR3.ends <- NULL
    }

    # --------------------------------------------------------------------------
    # UPSTREAM & DOWNSTREAM

    if ((!is.null(upstream)) | ((!is.null(downstream)))) {

      # Strand index
      idx.pos <- which(strand == "+")
      idx.neg <- which(strand == "-")

      # Chromosome length for every row; plus strand followed by minus strand
      if (!is.null(genome)) {
        chr.lengths <- unname(attr(genome, "length")[c(chrom[idx.pos],
                                                       chrom[idx.neg])])
      }
    }

    if (!is.null(upstream)) {

      # Upstream
      upstream.starts <- c(txStart[idx.pos] - upstream, txEnd[idx.neg] + 1)
      upstream.ends <- c(txStart[idx.pos] - 1, txEnd[idx.neg] + upstream)

      # Trim out-of-range genomic coordinate
      if (!is.null(genome)) {
        idx.upstream.ends.out <- which(upstream.ends > chr.lengths)

        upstream.starts[upstream.starts < 1] <- 1
        upstream.ends[idx.upstream.ends.out] <-
                               chr.lengths[idx.upstream.ends.out]
      }

    } else {
      upstream.starts <- upstream.ends <- NULL
    }

    if (!is.null(downstream)) {

      # Downstream
      downstream.starts <- c(txEnd[idx.pos] + 1, txStart[idx.neg] - downstream)
      downstream.ends <- c(txEnd[idx.pos] + downstream, txStart[idx.neg] - 1)

      # Trim out-of-range genomic coordinate
      if (!is.null(genome)) {
        idx.downstream.ends.out <- which(downstream.ends > chr.lengths)

        downstream.starts[downstream.starts < 1] <- 1
        downstream.ends[idx.downstream.ends.out] <-
                                         chr.lengths[idx.downstream.ends.out]
      }

    } else {
      downstream.starts <- downstream.ends <- NULL
    }

    # --------------------------------------------------------------------------
    # Combine and organise to a table

    if(!"exon" %in% element) exon.starts <- exon.ends <- NULL

    # Combine all coordinate above in a vector. The position arrangement is
    # important.
    start <- c(exon.starts, intron.starts, cds.starts, UTR5.starts, UTR3.starts,
               upstream.starts, downstream.starts)
    end   <- c(exon.ends, intron.ends, cds.ends, UTR5.ends, UTR3.ends,
               upstream.ends, downstream.ends)

    # Assign name, name2, chromosome, and strand to every coordinates above.
    # A function to help assign. The arrangement must be exactly the same like
    # above.
    repeat.feature <- function(feature){

      for.exons <- rep(feature, exonCount)
      for.introns <- if ("intron" %in% element) {
                       rep(feature, exonCount-1)
                     } else NULL
      for.cds <- if ("cds" %in% element) {
                   c(for.exons[idx.partly.cds.up & !idx.partly.one],
                     for.exons[idx.completely.cds],
                     for.exons[idx.partly.cds.down & !idx.partly.one],
                     for.exons[idx.partly.one])
                 } else NULL

      for.UTR5 <- if ("utr" %in% element) {
                    c(for.exons[idx.completely.UTR.up & idx.exon.pos],
                      for.exons[idx.partly.UTR.up & idx.exon.pos],
                      for.exons[idx.completely.UTR.down & idx.exon.neg],
                      for.exons[idx.partly.UTR.down & idx.exon.neg])
                  } else NULL

      for.UTR3 <- if ("utr" %in% element) {
                   c(for.exons[idx.partly.UTR.down & idx.exon.pos],
                     for.exons[idx.completely.UTR.down & idx.exon.pos],
                     for.exons[idx.partly.UTR.up & idx.exon.neg],
                     for.exons[idx.completely.UTR.up & idx.exon.neg])
                  } else NULL

      for.upstream <- if (!is.null(upstream)) {
                        c(feature[idx.pos], feature[idx.neg])
                      } else NULL
      for.downstream <- if (!is.null(downstream)) {
                          c(feature[idx.pos], feature[idx.neg])
                        } else NULL

      if (!"exon" %in% element) for.exons <- NULL

      return(c(for.exons, for.introns, for.cds, for.UTR5, for.UTR3,
               for.upstream, for.downstream))
    }

    name <- repeat.feature(name)
    if (type == "refseq") name2 <- repeat.feature(name2)
    chromosome <- repeat.feature(chrom)
    strand <- repeat.feature(strand)

    # Add type: "exon", "intron", "cds", "UTR5", "UTR3", "upstream",
    #           "downstream"
    element <- c(rep("exon", length(exon.starts)),
                 rep("intron", length(intron.starts)),
                 rep("CDS", length(cds.ends)),
                 rep("UTR5", length(UTR5.ends)),
                 rep("UTR3", length(UTR3.ends)),
                 rep("upstream", length(upstream.starts)),
                 rep("downstream", length(downstream.starts)))

    list(chromosome = chromosome, start = start, end = end, strand = strand,
         name = name, name2 = (if (type == "refseq") name2), element = element)
  }]

  # For noncoding transcripts, change 5'-UTR or 3'-UTR to UTR
  if (type == "refseq") {

    genomic.coordinate[ stri_sub_all(name, 2, 2)=="R" &
                       element %in% c("UTR5", "UTR3"), {

      element := "UTR"

    }, by = .(name, chromosome)]

  }

  # Convert back to 0-based indexing (-1 to every start coordinate)
  genome.annotation[, c("txStart", "cdsStart") := list(txStart-1, cdsStart-1)]

  if (nrow(genomic.coordinate) == 0) {
    message(paste("Elements [", element, "] are not found in the annotation",
                  "table."))
  }

  return(genomic.coordinate)
}
