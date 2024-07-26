#' Study k-mer composition across species.
#'
#' @param kmer.table A data.table of kmer table.
#' @param kmer.cutoff Percentage of extreme kmers to study. Default to 5.
#' @param genome.name UCSC genome name.
#' @param k K-mer size.
#' @param db Database used by UCSC to generate gene prediction: "refseq" or
#'     "gencode". Default is "refseq".
#' @param central.pattern K-mer's central patterns. Default is NULL.
#' @param output.dir A directory for the outputs.
#' @param fasta.path Path to a directory of user-provided genome FASTA files or 
#'  the destination to save the NCBI/UCSC downloaded reference genome files.
#' 
#' @return An output directory containing plots.
#'
#' @importFrom data.table fread fwrite setorder setnames rbindlist
#' @importFrom stringi stri_sub stri_count_regex stri_length stri_paste stri_sub_replace_all stri_replace_all_regex stri_split_fixed
#' @importFrom Biostrings reverseComplement
#' @importFrom graphics layout boxplot plot points legend arrows axis mtext barplot
#' @importFrom grDevices pdf
#' 
#' @export
STUDY_GENIC_ELEMENTS <- function(kmer.table, kmer.cutoff=5, k,
                                 genome.name="hg38", central.pattern=NULL,
                                 db="refseq",
                                 output.dir="study_genic_elements/",
                                 fasta.path) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  if (is.character(kmer.table))
    kmer.table <- fread(kmer.table, showProgress = FALSE)

  # Generate coordinate table for intergenic regions. --------------------------

  # Retrieve genePred table from UCSC.
  genepred <- getUCSCgenePredTable(genome.name, db)

  gene.name.col <- if(db == "refseq") "name2" else "geneName"

  # Filter genePred to include only chromosome references.
  genepred <- genepred[chrom %in% paste0("chr", c(1:22, "X", "Y"))]

  # Generate intergenic regions as a control.
  igr.dt <- generateIntergenicCoor(genepred, genome.name = genome.name, 
                                   fasta.path = fasta.path,
                                   igr.rel.pos = c(5000, 7500),
                                   igr.min.length = 150,
                                   return.coor.obj = FALSE)

  # Generate coordinate table for genic elements. ------------------------------

  # Only select complete CDS status
  genepred <- genepred[cdsStartStat == "cmpl" & cdsEndStat == "cmpl"]

  # Further restrict the gene to mRNA-coding only.
  genepred <- genepred[

    # Only include verified mRNA-coding transcript.
    if (db == "refseq") name %like% "NM_"
    else if (db == "gencode") transcriptClass == "coding"

  ]

  # Generate genic elements.
  gene.dt <- generateGenicElementCoor(genepred, element.names = "all",
                                      upstream = 2500, downstream = 1000,
                                      genome.name = genome.name,
                                      return.coor.obj = FALSE)

  # Remove exon as we are using UTR and CDS
  gene.dt <- gene.dt[element != "exon"]

  # Only take canonical gene construct i.e. must have both 5'-UTR and 3'-UTR
  #  1) 5'-UTR, CDS, intron, 3'-UTR
  #  2) 5'-UTR, CDS, 3'-UTR
  gene.dt <- gene.dt[, if(all(c("5'-UTR", "3'-UTR") %in% element)) .SD,
                     by = name]

  # Only include the longest transcript for alternative splice regions.
  # This is not chromosome-wise operation, so for similar genes in different
  # chromosome, only the longest one is selected.
  gene.dt <- gene.dt[, {

    sel.name <- .SD[, .(len = max(end) - min(start)),
                    by = name][which.max(len), name]

    .SD[name == sel.name][

      # This is for chrX and chrY
      if(length(unique(chromosome)) > 1) chromosome == chromosome[1] else TRUE]

  }, by = gene.name.col]

  # Combine tables of genic element and intergenic region.
  gene.dt <- rbind(gene.dt, igr.dt, fill = TRUE)

  # Consider one IGR region as a "gene"
  gene.dt[element == "IGR", (gene.name.col) := paste0("IGR_", 1:.N)]

  # Count k-mers ---------------------------------------------------------------
  setorder(kmer.table, -z)
  top.kmers <- kmer.table[1:round(kmer.cutoff / 100 * .N), kmer]
  bottom.kmers <- kmer.table[(.N - round(kmer.cutoff / 100 * .N)):.N, kmer]

  genome <- loadGenome(genome.name, fasta.style = "UCSC", fasta.path = fasta.path)

  gene.dt[element == "IGR", strand := "*"]
  setorder(gene.dt, chromosome, strand)

  flank <- (k - nchar(central.pattern)) / 2

  # Expand terminal
  gene.dt[, `:=`(start = start - flank, end = end + flank)]

  p <- progressr::progressor(steps = gene.dt[end - start + 1 >= k, 1,
    by = c("chromosome", "strand", gene.name.col, "element")]$V1 |> sum())

  t1 <- Sys.time()

  counts <- gene.dt[end - start + 1 >= k, {

    #message(paste(chromosome, strand, get(gene.name.col), element, "\n"))
    p(paste(chromosome, strand, get(gene.name.col), element))

    kmers <- buildRangedKmerTable(genome[chromosome], start, end, k,
                                  method = "chopping",
                                  chopping.method = "substring")

    # Temporarily name column N to pos_strand and create neg_strand column to
    # use countRevCompKmers.
    setnames(kmers, "N", "pos_strand")
    kmers[, neg_strand := 0]
    kmers <- countRevCompKmers(kmers)

    total.susc.pos <- stri_sub(genome[chromosome], start + flank, end -
      flank) |> stri_count_regex(central.pattern) |> sum()
    total.susc.neg <- stri_sub(genome[chromosome], start + flank, end -
      flank) |> stri_count_regex(reverseComplement(central.pattern)) |> sum()

    if (strand %in% c("+", "*")) {
      setnames(kmers, c("pos_strand", "neg_strand"), c("sense", "antisense"))
      total.susc.sense <- total.susc.pos
      total.susc.antisense <- total.susc.neg
    } else if (strand == "-") {
      setnames(kmers, c("neg_strand", "pos_strand"), c("sense", "antisense"))
      total.susc.sense <- total.susc.neg
      total.susc.antisense <- total.susc.pos
    }

    count <- kmers[
      stri_sub(kmer, flank + 1, flank + nchar(central.pattern)) ==
        central.pattern,

      .(X = c("sense", "antisense"),
        susceptible_kmer = c(total.susc.sense, total.susc.antisense),
        top_kmer = c(sum(sense[idx <- kmer %in% top.kmers]),
                     sum(antisense[idx])),
        bottom_kmer = c(sum(sense[idx <- kmer %in% bottom.kmers]),
                        sum(antisense[idx])))]
    count[, total_kmer := sum(end - start + 1 - k + 1)]

    for (col in names(count)[-1])
      set(count, j = col, value = as.double(count[[col]]))

    count

  }, by = c("chromosome", "strand", gene.name.col, "element")]

  time.taken <- Sys.time() - t1
  message(paste("Counting k-mers finished in ", round(time.taken, 2), " ",
      attr(time.taken, "units"), "\n", sep = ""))

  counts[, strand := NULL]
  setnames(counts, "X", "strand")

  fwrite(counts, paste0(output.dir, "/genic_elements_counts.csv"),
         showProgress = FALSE)

  # Plot -----------------------------------------------------------------------
  pdf(paste0(output.dir, "/plots.pdf"), paper = "a4",
      width = 8.27, height = 11.69)

  elm.cols <- c("khaki1", "cornflowerblue", "aquamarine4", "coral4",
                "dodgerblue4", "goldenrod3", "seashell4")
  names(elm.cols) <- c("upstream", "5'-UTR", "CDS", "intron", "3'-UTR",
                       "downstream", "IGR")

  counts[, `:=`(
    "susceptible k-mer content" = susceptible_kmer / total_kmer,
    "highly susceptible k-mer content" = top_kmer / susceptible_kmer
  )]
  # NaN because zero k-mer content. 0 / 0 = NaN
  counts[is.nan(`highly susceptible k-mer content`),
         "highly susceptible k-mer content" := 0]

  # Combine k-mer content in both strand
  counts.both <- counts[, .(
    susceptible_kmer = sum(susceptible_kmer),
    total_kmer = sum(total_kmer),
    top_kmer = sum(top_kmer)
  ), by = c(gene.name.col, "element")][, `:=`(
    "susceptible k-mer content" = susceptible_kmer / total_kmer,
    "highly susceptible k-mer content" = top_kmer / susceptible_kmer
  )]
  # NaN because zero k-mer content. 0 / 0 = NaN
  counts.both[is.nan(`highly susceptible k-mer content`),
         "highly susceptible k-mer content" := 0]

  # Boxplot
  layout(matrix(1:6, nrow = 3, byrow = FALSE))
  par(mar = c(7.1, 4.1, 4.1, 2.1))
  for (kmer.content in c("susceptible k-mer content",
                         "highly susceptible k-mer content")) {

    element.sorted <- names(elm.cols)
    # For both strands
    # element.sorted <- counts.both[, .(median = median(get(kmer.content),
    #                                                   na.rm = TRUE)),
    #                               by = element][order(median), c(element)]
    counts.both[, element := factor(element, levels = element.sorted)]
    boxplot(get(kmer.content) ~ element, data = counts.both,
            xaxt = "n", xlab = NA, col = elm.cols[element.sorted],
            boxcol = elm.cols[element.sorted],
            outline = FALSE, ylab = kmer.content, main = "both strands")
    axis(1, at = 1:7, labels = element.sorted, las = 2)
    title(xlab = "element", line = 6)

    # For separate strands
    # element.sorted <- counts[, .(median = median(get(kmer.content),
    #                                              na.rm = TRUE)),
    #                          by = element][order(median), c(element)]
    counts[, element := factor(element, levels = element.sorted)]

    counts[, {

      # element.sorted <- .SD[, .(median = median(get(kmer.content),
      #                                           na.rm = TRUE)),
      #                       by = element][order(median), c(element)]
      element <- factor(element, levels = element.sorted)
      boxplot(get(kmer.content) ~ element, outline = FALSE, xaxt = "n",
              main = paste0(strand, " strand"), ylab = kmer.content, xlab = NA,
              col = elm.cols[element.sorted], boxcol = elm.cols[element.sorted])
      axis(1, at = 1:7, labels = element.sorted, las = 2)
      title(xlab = "element", line = 6)

      NULL
    }, by = strand]
  }

  # Strand polarity plot
  # Order for calculating log2 fold change of strand i.e. strand polarity
  setorderv(counts, c("chromosome", gene.name.col, "element", "strand"))

  layout(matrix(1:2, nrow = 2, byrow = FALSE))
  par(mar = c(7.1, 5.1, 4.1, 2.1))
  for (kmer.content in c("susceptible k-mer content",
                         "highly susceptible k-mer content")) {

    # Calculate log2 sense / antisense
    fold.change <- counts[, .(
      element = element[strand == "sense"],
      log2_FC = log(get(kmer.content)[strand == "sense"] /
                      get(kmer.content)[strand == "antisense"],
                    base = 2))]

    # Convert NaN (0/0) as zero change.
    fold.change[is.na(log2_FC), log2_FC := 0]

    # Find the median/mean for building barplot. Ignore Inf fold change log(0)
    fc.med <- fold.change[is.finite(log2_FC),
      .(FC_median = mean(log2_FC), FC_sd = sd(log2_FC)), by = element]

    # Order from low to high value i.e. antisense to sense enrichment.
    setorder(fc.med, FC_median)
    fc.med[, element := factor(element, levels = names(elm.cols))]

    # Barplot
    axs.points <- barplot(FC_median ~ element, data = fc.med,
                          ylab = NA, xlab = NA, xaxt = "n",
                          main = paste0("Strand polarity of ", kmer.content),
                          col = elm.cols[fc.med[, levels(element)]],
                          border = elm.cols[fc.med[, levels(element)]])
    axis(1, at = axs.points, labels = names(elm.cols), las = 2, lwd = 0)
    title(xlab = "element", line = 6)
    title(ylab = expression("log"[2]*" ("*frac(sense, antisense)*")"),
          line = 2.5)
    # arrows(x0 = axs.points, x1 = axs.points,
    #        y0 = fc.med[match(levels(element), element), FC_median - FC_sd],
    #        y1 = fc.med[match(levels(element), element), FC_median + FC_sd],
    #        angle = 90, code = 3, length = 0.02, col = "gray")
  }

  dev.off()

}
