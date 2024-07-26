#' Study k-mer composition of selected COSMIC causal cancer genes across human
#'    populations worldwide.
#'
#' Simulation of human population is based on single nucleotide variantion.
#'
#' @param kmer.table A data.table of kmer table.
#' @param kmer.cutoff Percentage of extreme kmers to study. Default to 5.
#' @param genome.name UCSC genome name.
#' @param k K-mer size.
#' @param db Database used by UCSC to generate gene prediction: "refseq" or
#'     "gencode". Default is "refseq".
#' @param central.pattern K-mer's central patterns. Default is NULL.
#' @param population.size Size of population to simulate. Default is 1 million.
#' @param selected.genes Set of genes to study e.g. skin cancer genes.
#' @param output.dir A directory for the outputs. Default to
#'     study_across_populations.
#' @param add.to.existing.population Add counts to counts.csv? Default is
#'     FALSE.
#' @param population.snv.dt Population SNV table.
#' @param loop.chr Loop chromosome?. Default is TRUE. If FALSE, beware of a
#'     memory spike because of VCF content. VCF contains zero counts for every
#'     population. Input pre-computed trimmed-version population.snv.dt.
#' @param plot Boolean. Default is FALSE. If TRUE, will plot results.
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
STUDY_ACROSS_POPULATIONS <- function(kmer.table, kmer.cutoff=5, genome.name, k,
                                     db="refseq", central.pattern=NULL,
                                     population.size=1e6, selected.genes,
                                     add.to.existing.population=FALSE,
                                     output.dir="study_across_populations/",
                                     population.snv.dt=NULL, loop.chr=TRUE,
                                     plot=FALSE, fasta.path) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (is.character(kmer.table))
    kmer.table <- fread(kmer.table, showProgress = FALSE)

  if (is.character(population.snv.dt))
    population.snv.dt <- fread(population.snv.dt, showProgress = FALSE)

  gene.name.col <- if(db == "refseq") "name2" else "geneName"

  # Generate coordinate table for genic elements. ------------------------------

  # Retrieve genePred table from UCSC.
  genepred <- getUCSCgenePredTable(genome.name, db)

  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # Filter genePred to include only chromosome references.
  genepred <- genepred[paste0("chr", c(1:22, "X", "Y")) %in% chrom]

  # Only select selected genes
  genepred <- genepred[get(gene.name.col) %in% selected.genes]

  # Only select complete CDS status
  genepred <- genepred[cdsStartStat == "cmpl" & cdsEndStat == "cmpl"]

  # Prioritise NM_ over XM_ (This is only applicable to db = refseq)
  genepred <- genepred[, if(all(c("NM", "XM") %in% unique(pre <- stri_sub(name,
    1, 2)))) .SD[pre == "NM"] else .SD, by = eval(gene.name.col)]

  # Check if any selected genes is not in genePred
  missing.genes <- selected.genes[!selected.genes %in%
                                    genepred[, get(gene.name.col)]]
  if (length(missing.genes) > 0)
    warning("The following selected genes are not found in genePred: ",
            missing.genes)

  # Get genic elements
  # Column name is unique for every row, so the function is correct.
  gene.dt <- generateGenicElementCoor(genepred, element.names = "all",
                                      upstream = 2500, downstream = 1000,
                                      genome.name = genome.name,
                                      return.coor.obj = FALSE)

  # Remove exon as we are using UTR and CDS
  gene.dt <- gene.dt[element != "exon"]

  # Only take canonical gene construct i.e. must have both 5'-UTR and 3'-UTR
  #  1) 5'-UTR, CDS, intron, 3'-UTR
  #  2) 5'-UTR, CDS, 3'-UTR
  gene.dt <- gene.dt[, if (all(c("5'-UTR", "3'-UTR") %in% element)) .SD,
                     by = name]

  # Only include the longest transcript for alternative splice regions.
  # This is not chromosome-wise operation, so for similar genes in different
  # chromosome, only the longest one is selected.
  gene.dt <- gene.dt[, {

    sel.name <- .SD[, .(len = max(end) - min(start)),
                    by = name][which.max(len), name]

    .SD[name == sel.name]

  }, by = gene.name.col]

  # Check if any missing selected genes after filtering.
  missing.genes <- selected.genes[!selected.genes %in%
                                    gene.dt[, get(gene.name.col)]]
  if (length(missing.genes) > 0)
    warning("The following selected genes are removed after filtering step: ",
            missing.genes)

  fwrite(gene.dt, paste0(output.dir, "/genic_coordinates.csv"),
         showProgress = FALSE)

  # Expand the region to get k-meric context of the terminal base.
  flank <- (k - nchar(central.pattern)) / 2
  gene.dt[, `:=`(start = start - flank, end = end + flank)]

  # Identify extreme k-mers ----------------------------------------------------
  setorder(kmer.table, -z)
  top.kmers <- kmer.table[1:round(kmer.cutoff / 100 * .N), kmer]

  genome <- loadGenome(genome.name, fasta.style = "UCSC", fasta.path = fasta.path)

  pop.names <- c("afr", "ami", "amr", "asj", "eas", "fin", "nfe", "sas", "oth")

  element.names <- c("upstream", "5'-UTR", "CDS", "intron", "3'-UTR",
                     "downstream")

  # Delete these files if exist because we want to append
  if (is.null(population.snv.dt)) {
    unlink(paste0(output.dir, "/gnomad_SNV.vcf"))
    unlink(paste0(output.dir, "/population_SNV.csv"))
  }

  # Create table for k-mer counts in elements for each populations.
  counts <- as.data.table(matrix(0.0, nrow = population.size,
                                 ncol = length(pop.names) * 2 * 2 *
                                   length(element.names)))
  cols <- sapply(sapply(sapply(element.names, paste0,
                               c("_sense", "_antisense")),
                 paste0, "_", pop.names), paste0, "_", c("top", "susc"))
  setnames(counts, names(counts), cols)

  # Ignore strand for sorting because SNV info is on reference strand only.
  setorder(gene.dt, chromosome, start, end)

  # Start to operate element wise because consume a lot of memory!
  gene.dt[ifelse(population.size == 0 & !is.null(population.snv.dt), FALSE,
    TRUE), {

    if (.N == 0) return(NULL)

    # Merge regions to avoid VCF duplicate records; strand insensitive.
    merged.reg <- mergeCoordinate(.SD[, -"strand"])

    setkeyv(merged.reg, c(if(!loop.chr) "chromosome", "start", "end"))

    t1 <- Sys.time()

    if (is.null(population.snv.dt)) {

      # Get VCF in every elements - only variants pass the gnomAD QC pipeline.
      vcf <- merged.reg[, getGnomADvariants(chromosome, start, end,
        INFO.filter = c(paste0("AC_", pop.names), paste0("AN_", pop.names)))][
          FILTER == "PASS"]

      time.taken <- Sys.time() - t1
      message(paste("Finish parsing VCF"))
      message(paste("...", round(time.taken, digits = 2), " ",
          attr(time.taken, "units"), "\n", sep = ""))

      setorder(vcf, CHROM, POS)

      # Somehow set(vcf, j = "INFO", value = NULL) caused a segfault
      # writeVCF(vcf, paste0(output.dir, "/gnomad_SNV.vcf"), append = TRUE,
      #           tabix = FALSE)

      # The reason I want to writeVCF is because I want to explore, so these
      # filter will be excluded in population_SNV.csv
      # # Only select single base variant in REF and at least one in ALT
      # vcf <- vcf[stri_length(REF) == 1 &
      #              sapply(ALT, function(alt) any(stri_length(alt) == 1),
      #                     USE.NAMES = FALSE)]

      # Construct population frequency table with columns:
      # population, position, ID, chromosome, base, count
      pop.freq <- vcf[, {

        # Count
        alt.count <- unlist(.SD[, paste0("INFO.AC_", pop.names), with = FALSE],
                            use.names = FALSE)
        ref.count <- unlist(.SD[, paste0("INFO.AN_", pop.names), with = FALSE],
                            use.names = FALSE) -
          unlist(lapply(.SD[, paste0("INFO.AC_", pop.names), with = FALSE],
                        lapply, sum),
                 use.names = FALSE)
        count <- c(alt.count, ref.count)

        alt.lens <- lengths(ALT)

        ID <- c(rep(rep(ID, alt.lens), length(pop.names)),
                rep(ID, length(pop.names)))
        chromosome <- c(rep(rep(CHROM, alt.lens), length(pop.names)),
                        rep(CHROM, length(pop.names)))
        position <- c(rep(rep(POS, alt.lens), length(pop.names)),
                      rep(POS, length(pop.names)))
        population <- c(rep(pop.names, each = sum(alt.lens)),
                        rep(pop.names, each = .N))
        base <- c(rep(unlist(ALT, use.names = FALSE), length(pop.names)),
                  rep(REF, length(pop.names)))
        REF <- c(rep(rep(REF, alt.lens), length(pop.names)),
                 rep(REF, length(pop.names)))

        .(ID = ID, chromosome = chromosome, position = position,
          population = population, ref_base = REF, base = base, count = count)
      }]

      rm(vcf)

      setorder(pop.freq, chromosome, position)
      fwrite(pop.freq, paste0(output.dir, "/population_SNV.csv"),
             append = TRUE, showProgress = FALSE)
    } else {

      pop.freq <- population.snv.dt[chromosome %in% unique(chromosome)]
      setorder(pop.freq, chromosome, position)
      pop.freq[, position2 := position]
      pop.freq <- foverlaps(pop.freq, merged.reg, by.x = c(if(!loop.chr)
        "chromosome", "position", "position2"), nomatch = NULL)
      pop.freq[, c("start", "end", "position2",
                   names(pop.freq)[grepl("^i\\.", names(pop.freq))]) := NULL]
    }

    # Remove base duplicates.
    # Observation: Duplicate position; one with rsID and one without.
    # Remove duplicates (count prioritised) because it will throw off the
    # probabilities. It seems the major allele raw number are about the same?
    # Hmmmm.... Are they from the same study? Need more investigation from the
    # saved population_SNV.csv.
    setorder(pop.freq, population, chromosome, position, base, -count)
    pop.freq <- pop.freq[

      # Remove non-SNV, if any.
      stri_length(base) == 1 &

      # Remove zero count and base duplicates with lower count.
      count != 0,

      .(ref_base = ref_base[1], count = count[1], ID = ID[1]),
      by = .(population, chromosome, position, base)][

        # Only take population with multiple alleles or if one allele, it must
        # be different from the reference base.
        , if (.N > 1 || stri_sub(ref_base, 1, 1) != base) .SD,
        by = .(population, chromosome, position)]

    # Append merged region sequences
    merged.seq <- merged.reg[, .(stri_sub(genome[chromosome], start, end) |>
      stri_paste(collapse = "N")), by = eval(if(!loop.chr) "chromosome")]$V1 |>
      stri_paste(collapse = "N")

    # Resolve SNV relative positions; strand insensitive
    merged.reg[, group := 1:.N - 1]
    merged.reg[, `:=`(rel_end = end - start + 1)]
    merged.reg[, shift := cumsum(shift(rel_end, fill = 0)) + group]
    merged.reg[, rel_end := NULL]
    pop.freq[, position2 := position]
    pop.freq <- foverlaps(pop.freq, merged.reg, by.x = c(if(!loop.chr)
      "chromosome", "position", "position2"))
    pop.freq[, position := position - start + 1 + shift]
    pop.freq[, c("start", "end", "group", "shift", "position2",
                 names(pop.freq)[grepl("^i\\.", names(pop.freq))]) := NULL]

    # Resolve element coordinates
    rel.coors <- foverlaps(.SD, merged.reg)
    rel.coors[, `:=`(start = i.start - start + 1 + shift,
                     end = i.end - start + 1 + shift)]
    rel.coors[, c("group", "shift",
                  names(rel.coors)[grepl("^i\\.", names(rel.coors))]) := NULL]

    p <- progressr::progressor(steps = population.size * length(pop.names))

    # Simulation of population individuals.
    pop.freq[if(population.size == 0) FALSE else TRUE, {

      t1 <- Sys.time()

      counts.pop <- future.apply::future_replicate(n = population.size, {

        p(paste0("simulation: ", if(loop.chr) chromosome, population))

        # t1.a <- Sys.time()

        # Vectorized sampling
        snv <- cbind(.SD, prob = runif(.N) * count) |> setorder(position, -prob)
        snv <- snv[, .(base = base[1]), by = position]

        # t1.b <- Sys.time()
        # time.taken <- t1.b - t1.a
        # message(paste("Sampling bases"))
        # message(paste("...", round(time.taken, digits = 2), " ",
        #     attr(time.taken, "units"), "\n", sep = ""))

        mut.merged.seq <- stri_sub_replace_all(merged.seq,
                                               from = snv$position,
                                               length = stri_length(snv$base),
                                               value = snv$base)

        # t1.c <- Sys.time()
        # time.taken <- t1.c - t1.b
        # message(paste("Mutating genome"))
        # message(paste("...", round(time.taken, digits = 2), " ",
        #     attr(time.taken, "units"), "\n", sep = ""))

        # Count susceptible k-mers in every element for each strand.
        i.cnts <- data.table()
        rel.coors[, {

          # t1.b <- Sys.time()

          total.susc.pos <- stri_sub(mut.merged.seq, start + flank, end -
            flank) |> stri_count_regex(central.pattern) |> sum()
          total.susc.neg <- stri_sub(mut.merged.seq, start + flank, end -
            flank) |> stri_count_regex(reverseComplement(central.pattern)) |>
            sum()

          kmers <- buildRangedKmerTable(mut.merged.seq, start, end, k,
                                        method = "chopping",
                                        chopping.method = "Biostrings",
                                        remove.N = TRUE)

          # Count on the opposite strand
          setnames(kmers, "N", "pos_strand")
          kmers[, neg_strand := 0]
          kmers <- countRevCompKmers(kmers)

          if (strand %in% c("+", "*")) {
            setnames(kmers, c("pos_strand", "neg_strand"),
                     c("sense", "antisense"))
            total.susc.sense <- total.susc.pos
            total.susc.antisense <- total.susc.neg
          } else if (strand == "-") {
            setnames(kmers, c("neg_strand", "pos_strand"),
                     c("sense", "antisense"))
            total.susc.sense <- total.susc.neg
            total.susc.antisense <- total.susc.pos
          }

          # Set all susceptible k-mers
          set(i.cnts, j = paste0(element, "_sense_", population, "_susc"),
              value = as.numeric(total.susc.sense))
          set(i.cnts, j = paste0(element, "_antisense_", population, "_susc"),
              value = as.numeric(total.susc.antisense))

          # Count relevant k-mers
          kmers[kmer %in% top.kmers, {

            set(i.cnts,
                j = paste0(element, "_sense_", population, "_top"),
                value = as.numeric(sum(sense)))
            set(i.cnts,
                j = paste0(element, "_antisense_", population, "_top"),
                value = as.numeric(sum(antisense)))

            NULL
          }]

          # t1.c <- Sys.time()
          # time.taken <- t1.c - t1.b
          # message(paste("Counting k-mers in", element, "of", strand, "strand"))
          # message(paste("...", round(time.taken, digits = 2), " ",
          #     attr(time.taken, "units"), "\n", sep = ""))

        }, by = .(element, strand)]

        # t2.a <- Sys.time()
        # time.taken <- t2.a - t1.a
        # message(paste("Individual from", population, "population, zooming into",
        #     "all elements of selected genes", if(loop.chr) c("in", chromosome)))
        # message(paste("...", round(time.taken, digits = 2), " ",
        #     attr(time.taken, "units"), "\n", sep = ""))

        i.cnts
      }, future.globals = c(".N", ".SD", "count", "merged.seq", "rel.coors",
                            "central.pattern", "top.kmers", "flank", "element",
                            "population", "loop.chr", "chromosome", "strand",
                            "p"),
      simplify = FALSE, future.seed = TRUE) |> rbindlist()

      counts[, names(counts.pop) := .SD + counts.pop,
             .SDcols = names(counts.pop)]

      t2 <- Sys.time()
      time.taken <- t2 - t1
      message(paste("Simulation of", population.size, "individuals of", population,
          "population in all elements", if(loop.chr) c("of", chromosome)))
      message(paste("...", round(time.taken, digits = 2), " ",
          attr(time.taken, "units"), "\n", sep = ""))

      NULL
    }, by = eval(c(if(loop.chr) "chromosome", "population"))]

    NULL
  }, by = eval(if(loop.chr) "chromosome")]

  # Save total k-mers in a single strand
  gene.dt[, {
    counts[, paste0(element, "_total_kmers") := sum(end - start + 1 - k + 1)]
    NULL
  }, by = element]

  fwrite(counts, paste0(output.dir, "/counts.csv"),
         append = add.to.existing.population,
         showProgress = FALSE)

  if (!plot) return(counts)

  # Plot -----------------------------------------------------------------------
  pdf(paste0(output.dir, "/plots.pdf"), paper = "a4r",
      width = 11.69, height = 8.27)

  pop.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#FFFF33", "#A65628", "#F781BF", "#999999") |> addAlphaCol(0.9)
  names(pop.cols) <- pop.names

  # Load population frequency
  if (is.data.table(population.snv.dt)) {
    pop.freq <- population.snv.dt
  } else {
    pop.freq <- fread(paste0(output.dir, "/population_SNV.csv"),
                    showProgress = FALSE)
  }
  setorder(pop.freq, population, chromosome, position, base, -count)
  pop.freq <- pop.freq[

    # Remove non-SNV, if any.
    stri_length(base) == 1 &

    # Remove zero count and base duplicates with lower count.
    count != 0,

    .(ref_base = ref_base[1], count = count[1], ID = ID[1]),
    by = .(population, chromosome, position, base)][

      # Only take population with multiple alleles or if one allele, it must
      # be different from the reference base.
      , if (.N > 1 || stri_sub(ref_base, 1, 1) != base) .SD,
      by = .(population, chromosome, position)]

  # SNV position distribution
  pop.freq[, .(N = length(unique(position))), by = population][order(N), {

    par(mar = c(5.1, 4.1, 4.1, 4.1))
    barplot(N, names = population, col = pop.cols[population],
            xlab = "population", ylab = "SNV-position count", border = NA,
            main = "SNV distribution")
    par(new = TRUE)
    barplot(N / sum(N) * 100, col = NA, border = NA, xaxt = "n", yaxt = "n")
    axis(side = 4)
    mtext("percentage", side = 4, line = 3)

  }]

  # Venn diagram of population SNV position

  # Plot count distribution
  counts <- fread(paste0(output.dir, "/counts.csv"), showProgress = FALSE)

  # Add both strand count
  sense.cols <- names(counts)[grep("_sense_", names(counts))]
  counts[, sub("_sense_", "_both_", sense.cols) :=
           .SD[, sense.cols, with = FALSE] +
           .SD[, sub("_sense_", "_antisense_", sense.cols), with = FALSE]]

  # Plot total counts
  counts[, {

    # Plot layout 3 x 3 - 6 boxplots and 3 scatter plots
    lay.mat <- matrix(1:9,
                      nrow = 3, byrow = FALSE)
    layout(lay.mat)

    tot.cnt <- rowSums(.SD[, grep(paste0("_total_kmers"), names(.SD)),
                              with = FALSE]) |> c() |> unique()
    if (length(tot.cnt) > 1) stop("Total all k-mers should be the same!")

    for (strand.type in c("both", "sense", "antisense")) {
      top.cnts <- sapply(pop.names, function(pop.name) {
        rowSums(.SD[, grep(paste0(strand.type, "_", pop.name, "_top"),
          names(.SD)), with = FALSE])
      })
      susc.cnts <- sapply(pop.names, function(pop.name) {
        rowSums(.SD[, grep(paste0(strand.type, "_", pop.name, "_susc"),
          names(.SD)), with = FALSE])
      })

      # Convert to proportion
      top.cnts <- top.cnts / susc.cnts
      top.cnts[is.nan(top.cnts)] <- 0
      susc.cnts <- susc.cnts / (tot.cnt * ifelse(strand.type == "both", 2, 1))

      # Boxplot
      boxplot(susc.cnts, col = pop.cols[colnames(susc.cnts)],
              boxcol = pop.cols[colnames(susc.cnts)], outline = FALSE,
              main = paste0("Susceptible k-mers in ", strand.type, " strand",
                            if(strand.type == "both") "s"),
              xlab = "k-mer content")
      boxplot(top.cnts, col = pop.cols[colnames(top.cnts)],
              boxcol = pop.cols[colnames(top.cnts)], outline = FALSE,
              main = paste0("Highly susceptible k-mers in ", strand.type,
                            " strand", if(strand.type == "both") "s"),
              xlab = "k-mer content")

      # Scatter plot
      top.means <- mean(top.cnts)
      susc.means <- mean(susc.cnts)
      top.sdevs <- sd(top.cnts)
      susc.sdevs <- sd(susc.cnts)

      # Initialise empty plot.
      plot(x = susc.means, y = top.means, type = "n",
           xlim = range(susc.means - susc.sdevs, susc.means + susc.sdevs),
           ylim = range(top.means - top.sdevs, top.means + top.sdevs),
           xlab = "Susceptible k-mer content",
           ylab = "Highly susceptible k-mer content")

      # Standard deviation in susceptible k-mer content
      arrows(x0 = susc.means - susc.sdevs, y0 = top.means,
             x1 = susc.means + susc.sdevs, y1 = top.means,
             angle = 90, code = 3, length = 0.02, col = "gray")

      # Standard deviation in highly-susceptible k-mer content
      arrows(x0 = susc.means, y0 = top.means - top.sdevs,
             x1 = susc.means, y1 = top.means + top.sdevs,
             angle = 90, code = 3, length = 0.02, col = "gray")

      points(x = susc.means, y = top.means, pch = 16, col = pop.cols[pop.names],
             cex = 1.2)
    }
  }]

  # Plots by elements
  # Convert count to proportion
  # top count / susc count
  top.cols <- names(counts)[grep("_top$", names(counts))]
  counts[, c(top.cols) := .SD[, top.cols, with = FALSE ] /
           .SD[, sub("_top$", "_susc", top.cols), with = FALSE]]

  # susc count / total k-mers
  elements <- unique(stri_split_fixed(names(counts), "_", simplify = TRUE)[, 1])
  for (element in elements) {
    count.cols <- names(counts)[grep(paste0(element, "_(sense|antisense)_susc"),
                                     names(counts))]
    counts[, c(count.cols) := .SD[, count.cols, with = FALSE] /
      .SD[, paste0(element, "_total_kmers"), with = FALSE]]
    count.cols <- names(counts)[grep(paste0(element, "_both_susc"),
                                     names(counts))]
    counts[, c(count.cols) := .SD[, count.cols, with = FALSE] /
      (.SD[, paste0(element, "_total_kmers"), with = FALSE] * 2)]
  }

  # NaN because zero k-mer content. 0 / 0 = NaN. Convert to zero.
  for (col in names(counts)) {
    set(counts, i = which(is.nan(counts[[col]])), j = col, value = 0)
  }

  means <- sapply(counts, mean)
  sdevs <- sapply(counts, sd)

  # Plot layout.
  # Panel 1-6 are scattered plots.
  # Panel 7 is title.
  # Panel 8 is gene construct drawing.
  # Panel 9 is ylab
  # Panel 10 is xlab
  lay.mat <- matrix(c(0, rep(7, 3),
                      0, 1:3,
                      0, rep(8, 3),
                      9, 4:6,
                      0, 10, 0, 0),
                    nrow = 5, byrow = TRUE)
  layout(lay.mat, heights = c(2, 5, 3, 5, 0.875),
         widths = c(0.5, 4, 4, 4))

  par(mgp = c(3, 0.7, 0)) # Move tick labels close to the tick.

  # Order is important.
  element.names <- c("upstream", "CDS", "3'-UTR", "5'-UTR", "intron",
                     "downstream")

  for (strand.type in c("both", "sense", "antisense")) {

    # Panel 1-6 scatter plots
    par(mar = c(1, 1, 1, 1))

    for (element in element.names) {

      top.props <- means[grep(paste0(element, "_", strand.type, "_.*_top$"),
                              names(means))]
      susc.props <- means[sub("_top$", "_susc", names(top.props))]

      top.sdevs <- sdevs[names(top.props)]
      susc.sdevs <- sdevs[names(susc.props)]

      # Initialise empty plot.
      plot(x = susc.props, y = top.props, type = "n",
           xlim = range(susc.props - susc.sdevs, susc.props + susc.sdevs),
           ylim = range(top.props - top.sdevs, top.props + top.sdevs),
           axes = FALSE, ann = FALSE)
      axis(1, col = NA, col.ticks = "azure3")
      axis(2, col = NA, col.ticks = "azure3")
      box(col = "azure3")
      mtext(element, at = par("usr")[c(1, 4)], adj = 0, cex = 0.6,
            col = "azure3")

      # Standard deviation in susceptible k-mer content
      arrows(x0 = susc.props - susc.sdevs, y0 = top.props,
             x1 = susc.props + susc.sdevs, y1 = top.props,
             angle = 90, code = 3, length = 0.02, col = "gray")

      # Standard deviation in highly-susceptible k-mer content
      arrows(x0 = susc.props, y0 = top.props - top.sdevs,
             x1 = susc.props, y1 = top.props + top.sdevs,
             angle = 90, code = 3, length = 0.02, col = "gray")

      pop.order <- stri_extract_first_regex(names(top.props),
                                            paste0("_(", paste(pop.names,
                                                               collapse = "|"),
                                                   ")_")) |>
        stri_replace_all_regex("_", "")
      points(x = susc.props, y = top.props, pch = 16, col = pop.cols[pop.order],
             cex = 1.2)

    }

    # Panel 7 - title
    par(mar = c(0, 1, 0, 1))
    plot(1, type = "n", axes = FALSE, ann = FALSE)
    mtext(paste0("K-mer content in ", strand.type, " strand"), side = 3,
          line = -3, font = 2)
    legend(x = par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2,
           y = par("usr")[3] + (par("usr")[4] - par("usr")[3]) / 2,
           legend = names(pop.cols), col = pop.cols, bty = "n", horiz = TRUE,
           pch = 16, xjust = 0.5, cex = 1.2)

    # Panel 8 - Gene construct drawing
    par(mar = c(0, 1, 2, 1))
    plot(NULL, xlim = c(0, 11), ylim = c(0, 10), xaxs = "i", yaxs = "i",
         axes = FALSE, ann = FALSE)

    # Upstream
    lines(c(0, 3), c(5, 5))
    lines(c(0, 0), c(4.7, 5.3))

    # Symbol for the exon beginning
    lines(c(3, 3), c(6, 7), lwd = 0.7, col = "dimgrey")
    lines(c(3, 3.25), c(7, 7), lwd = 0.7, col = "dimgrey")
    polygon(c(3.25, 3.25, 3.3), c(6.7, 7.3, 7), lwd = 0.7, col = "dimgrey",
            border = NA)

    # 5'-UTR
    rect(3, 4, 4, 6, border = NA, col = "grey")

    # CDS
    rect(4, 4, 5, 6, border = NA, col = "seagreen")

    # Intron
    lines(c(5, 6), c(5, 7))
    lines(c(6, 7), c(7, 5))

    # CDS
    rect(7, 4, 8, 6, border = NA, col = "seagreen")

    # 3'-UTR
    rect(8, 4, 9, 6, border = NA, col = "grey")

    # Symbol for the exon ending
    points(9, 5, pch = 16)

    # Downstream
    lines(c(9, 11), c(5,5))
    lines(c(11, 11), c(4.7, 5.3))

    # Arrows pointing to plots
    arrows(1.5, 7.5, 1.5, 9.5, length = 0.05)
    arrows(3.5, 3, 3, 1, length = 0.05)
    arrows(4.5, 7.5, 4.5, 9.5, length = 0.05)
    arrows(6, 3, 6, 1, length = 0.05)
    arrows(8.5, 7.5, 8.5, 9.5, length = 0.05)
    arrows(10, 3, 10, 1, length = 0.05)

    # Panel 9 - ylab
    par(mar = c(1, 1, 1, 0))
    plot(1, type = "n", axes = FALSE, ann = FALSE)
    mtext("highly-susceptible k-mer content", side = 4, line = -2.5, cex = 0.8)

    # Panel 10 - xlab
    par(mar = c(1, 1, 0, 1))
    plot(1, type = "n", axes = FALSE, ann = FALSE)
    mtext("susceptible k-mer content", side = 3, line = -2.5, cex = 0.8)

  }

  dev.off()

  # Idea: 3 heatmap of cross-pair p-values: sense, antisense, both

}