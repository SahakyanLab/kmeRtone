bootstrappingGenes <- function(transcript.susceptibility, cancer.gene.census,
                               tumour.keyword, tumour.exact.keyword,
                               seed.number = 1, replication.number,
                               subsample.number = NULL, output, ncpu = 1) {
  # Dependencies
  # Packages: data.table, R.utils

  # DEV - WIP
  # To use nice coloring

  # A function to make to calculate density based on the same x. As a result,
  # all density calculated will have the same x points on a plot.
  density.ref <- function(vec, x.ref) {
    density(vec, n = length(x.ref),
               from = min(x.ref),
               to = max(x.ref))
  }

  # A function to perform bootstrap.
  # 1. Sample the data for many events
  # 2. Calculate the density using density.ref function
  # 3. Because x is the same, only output the Y value (density) in matrix
  #                               replicate 1, replicate 2, 3, ...
  #    density at kmer fraction 1
  #    density at kmer fraction 2
  #                           ...
  calculate.bootstrap.density <- function(content, main.density) {

    operation <- function() {
      bootstrap.density <-
        replicate(replication.number,
                density.ref(vec = sample(content, subsample.number),
                            x.ref = main.density$kmer_fraction)$y)

      return(bootstrap.density)
    }

    if (ncpu == 1) {

      bootstrap.density <- operation()

    } else {

      replication.numbers <- distributeChunk(replication.number, ncpu)$size
      seed.numbers <- rep(replication.number, ncpu) + 1:ncpu

      to.export.variable <- c("seed.numbers", "subsample.number", "content",
                              "main.density")
      to.export.function <- c("density.ref", "operation")

      bootstrap.density <-
        foreach(replication.number = replication.numbers,
                seed.number = seed.numbers,
                .combine = "cbind", .packages = c("data.table"),
                .export = c(to.export.variable, to.export.function),
                .noexport = ls()[!ls() %in% to.export.variable]) %dopar% {

          set.seed(seed.number)
          bootstrap.density <- operation()
          return(bootstrap.density)
      }
    }
    return(bootstrap.density)
  }

  # A function to calculate -log(p-value) named as sig score in each kmer
  # fraction population. The output is a vector of sig score.
  calculate.p.sig.score <- function(the.density, bootstrap.density) {

    side <- ifelse(the.density$target_cancers > rowMeans(bootstrap.density),
                   "greater", "less")

    p.value <- mapply(wilcox.test,

      # parameters for wilcox.test function
      x = the.density$target_cancers,
      y = split(bootstrap.density, 1:nrow(bootstrap.density)),
      alternative = side,
      SIMPLIFY = FALSE
    )

    p.value <- sapply(p.value, function(i) i$p.value)

    p.sig.score <- -log(p.value)
    p.sig.score[side == "less"] <- p.sig.score[side == "less"] * -1

    return(p.sig.score)
  }

  # Shade significant area: enhance and suppress
  shade.sig <- function(sig.points, color, the.density) {

    # A matrix with two columns: from, to
    sig.areas <- seqToIntervals(sig.points, color)

    rect(xleft = the.density$kmer_fraction[sig.areas[, "from"]],
         xright = the.density$kmer_fraction[sig.areas[, "to"]],
         ybottom = 0, ytop = par("usr")[4],
         col = color,
         border = NA)

    # To draw or plot only in this region. The rest will be clipped.
    # This is to ensure abline is truncated and start from y = zero
    # DEV - x1 should be zero but need to trim the kmer fraction outliers
    #       first
    #clip(x1 = par("usr")[1], x2 = par("usr")[2], y1 = 0, y2 = par("usr")[4])

    # Draw only vertical border line for the rectangle
    #abline(v = the.density$kmer_fraction[as.vector(sig.areas)],
    #       col = color, lwd = 0.7)

  }

  # A plot function
  plot.me <- function(the.density, bootstrap.density, p.sig.score,
                      content.label, element, strand, score, strand.label) {
    # label   "base content", "high risk kmers", "low risk kmers"

    par(fig = c(0, 1, 0, 0.8))

    # Initiate plot
    plot(x = range(the.density$kmer_fraction),
         y = range(bootstrap.density, the.density[, -1]),
         type = "n", # Don't display anything
         xlab = "kmer fraction",
         ylab = "density")

    # Bootstrap of all genes
    matlines(the.density$kmer_fraction, bootstrap.density,
             col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.05),
             lty = 1)
    # All genes
    lines(x = the.density$kmer_fraction, y = the.density$all_genes,
         ylim = range(the.density[, -1], bootstrap.density),
         col =  c1 <- "#0066807f",
         lty = 1,
         lwd = 2
    )
    # Target cancer genes
    lines(the.density$kmer_fraction,
          the.density$target_cancers,
          col = c2 <- "#8000337f",
          lty = 1,
          lwd = 2)

    # Color background of p.sig.score
    sig.enhanced.points <- which(p.sig.score > -log(0.05))
    sig.suppressed.points <- which(p.sig.score < log(0.05))

    if (length(sig.enhanced.points) > 0)
      shade.sig(sig.enhanced.points, c3 <- "#a02c2c4d", the.density)
    if (length(sig.suppressed.points) > 0)
      shade.sig(sig.suppressed.points, c4 <- "#71c8374d", the.density)

    legend("topright",
           legend = c("all genes", "targeted cancer genes"),
           fill = c(c1, c2), cex = 0.8, box.lty = 0)

    # plot p.sig.score i.e. -log(p-value) on top
    par(fig = c(0, 1, 0.55, 1), new = TRUE)
    plot(the.density$kmer_fraction, p.sig.score, type = "l", xaxt='n',
         xlab = "", ylab = "significant score")

    # shade sig area on top sig.score graph
    if (length(sig.enhanced.points) > 0)
      rect(xleft = par("usr")[1], ybottom = -log(0.05),
           xright = par("usr")[2], ytop = par("usr")[4],
           col = c3, border = NA)
    if (length(sig.suppressed.points) > 0)
      rect(xleft = par("usr")[1], ybottom = par("usr")[3],
           xright = par("usr")[2], ytop = log(0.05),
           col = c4, border = NA)

    title(main = paste("Distribution of", content.label, "in", element, "of",
                      c(if (length(strand) > 1) "both" else strand),
                      "strand"))
  }

  # ----------------------------------------------------------------------------
  # MAIN code

  set.seed(seed.number)

  # Resolve regex keywords
  tumour.keyword <- keyword2regex(tumour.keyword, tumour.exact.keyword)

  # Change column name to match genome.annotation table for table merging
  setnames(cancer.gene.census, "Gene Symbol", "name2")

  # Add column "Tumour Types(Somatic)" and "Tumour Types(Germline)" to
  # transcript.susceptibility table
  transcript.susceptibility <-
    merge.data.table(transcript.susceptibility,
                     cancer.gene.census[, c("name2", "Tumour Types(Somatic)",
                                            "Tumour Types(Germline)")],
                     by = "name2", all.x = TRUE)

  # Default value for subsample.number if not supplied
  if (is.null(subsample.number)) {
    subsample.number <-
      transcript.susceptibility[grepl(tumour.keyword, `Tumour Types(Somatic)`) |
                                grepl(tumour.keyword, `Tumour Types(Germline)`),
                              length(unique(name2))]
  }

  # This command is specific to NCBIrefSeq table. I'm not sure about GENCODE
  setkey(transcript.susceptibility, strand, element, name, name2)

  # Scores in the table e.g. z, z_TT, z_AA, etc.
  scores <- names(transcript.susceptibility)[
                          grepl("risk|susceptible",
                                names(transcript.susceptibility))]
  scores <- unique(gsub("susceptible_kmers_|high_risk_kmers_|low_risk_kmers_",
                        "", scores))

  for (score in scores) {

    pdf(paste0(output, "/density_bootstrap_", score, ".pdf"))

    for (both.strand in c(TRUE, FALSE)) {

      transcript.susceptibility[element != "IGR", {

        base.content <- eval(parse(text = paste0("susceptible_kmers_", score))) /
                        total_kmers

        high.risk.content <- eval(parse(text = paste0("high_risk_kmers_",
                                                      score))) /
                             eval(parse(text = paste0("susceptible_kmers_",
                                                      score)))

        low.risk.content <- eval(parse(text = paste0("low_risk_kmers_",
                                                      score))) /
                            eval(parse(text = paste0("susceptible_kmers_",
                                                      score)))

        # strand.label for pdf name
        strand.label <- strand

        if (both.strand) {

          idx.sense <- which(strand == "sense")
          idx.antisense <- which(strand == "antisense")

          base.content <- (base.content[idx.sense] +
                           base.content[idx.antisense]) / 2

          high.risk.content <- (high.risk.content[idx.sense] +
                                high.risk.content[idx.antisense]) / 2

          low.risk.content <- (low.risk.content[idx.sense] +
                               low.risk.content[idx.antisense]) / 2
          strand.label <- "both_strand"
        }

        # Total number of data. This is purely for both strand dataset because
        # we reduce the dataset by half in the if condition above. It won't affect
        # the single strand operation even though it is applied globally.
        n <- length(base.content)

        # Index of non-NaN number for high or low risk.
        # NaN as a result of 0/0 where there is no susceptible kmers at all.
        # Any mathematical operation involving NaN results in NaN.
        idx.hrisk.finite <- !is.nan(high.risk.content)[1:n]
        idx.lrisk.finite <- !is.nan(low.risk.content)[1:n]

        # Take index of all cancers and target cancer
        idx.all.cancer <- ((!is.na(`Tumour Types(Somatic)`)) |
                           (!is.na(`Tumour Types(Germline)`)))[1:n]
        idx.target.cancer <- (grepl(tumour.keyword, `Tumour Types(Somatic)`) |
                             grepl(tumour.keyword, `Tumour Types(Germline)`))[1:n]

        # Calculate density

        base.content.density <- data.table(

          kmer_fraction = {ref <- density(base.content); ref$x},
          all_genes = ref$y,
          all_cancers = density.ref(base.content[idx.all.cancer], ref$x)$y,
          target_cancers = density.ref(base.content[idx.target.cancer], ref$x)$y
        )

        high.risk.density <- data.table(

          kmer_fraction = {ref <- density(high.risk.content[idx.hrisk.finite])
                           ref$x},
          all_genes = ref$y,
          all_cancers = density.ref(high.risk.content[idx.all.cancer &
                                      idx.hrisk.finite], ref$x)$y,
          target_cancers = density.ref(high.risk.content[idx.target.cancer &
                                         idx.hrisk.finite], ref$x)$y
        )

        low.risk.density <- data.table(

          kmer_fraction = {ref <- density(low.risk.content[idx.lrisk.finite])
                           ref$x},
          all_genes = ref$y,
          all_cancers = density.ref(low.risk.content[idx.all.cancer &
                                      idx.lrisk.finite], ref$x)$y,
          target_cancers = density.ref(low.risk.content[idx.target.cancer &
                                         idx.lrisk.finite], ref$x)$y
        )

        # Sanity check: Check x values to be the same for every replicate. [KIV]
        # p/s: Claudia did calculate non-skin cancer genes and plot the graph.
        #      I will revisit this later [DEV]. Note that I need to substract
        #      skin cancer genes from the all cancer genes.

        bootstrap.base.content.density <-
          calculate.bootstrap.density(base.content, base.content.density)

        bootstrap.hrisk.density <-
          calculate.bootstrap.density(high.risk.content[idx.hrisk.finite],
                                      high.risk.density)

        bootstrap.lrisk.density <-
          calculate.bootstrap.density(low.risk.content[idx.lrisk.finite],
                                      low.risk.density)

        # DEV
        # Get optimum range i.e. remove density outliers
        temp.range <- boxplot.stats(base.content)$stats[c(1, 5)]
        temp.idx <- base.content.density[kmer_fraction >= temp.range[1] &
                                         kmer_fraction <= temp.range[2],
                                       which = TRUE]
        base.content.density <- base.content.density[temp.idx]
        bootstrap.base.content.density <-
          bootstrap.base.content.density[temp.idx, ]

        temp.range <-
          boxplot.stats(high.risk.content[idx.hrisk.finite])$stats[c(1, 5)]
        temp.idx <- high.risk.density[kmer_fraction >= temp.range[1] &
                                      kmer_fraction <= temp.range[2],
                                    which = TRUE]
        high.risk.density <- high.risk.density[temp.idx]
        bootstrap.hrisk.density <- bootstrap.hrisk.density[temp.idx, ]

        temp.range <-
          boxplot.stats(low.risk.content[idx.lrisk.finite])$stats[c(1, 5)]
        temp.idx <- low.risk.density[kmer_fraction >= temp.range[1] &
                                     kmer_fraction <= temp.range[2],
                                   which = TRUE]
        low.risk.density <- low.risk.density[temp.idx]
        bootstrap.lrisk.density <- bootstrap.lrisk.density[temp.idx, ]

        # Calculate p-value for each kmer fraction population.
        p.base.content.target.cancer <-
          calculate.p.sig.score(base.content.density,
                                bootstrap.base.content.density)

        p.hrisk.target.cancer <-
          calculate.p.sig.score(high.risk.density, bootstrap.hrisk.density)

        p.lrisk.target.cancer <-
         calculate.p.sig.score(low.risk.density, bootstrap.lrisk.density)

        # Plot

        plot.me(base.content.density, bootstrap.base.content.density,
                p.base.content.target.cancer, "base content", element, strand,
                score, strand.label)
        plot.me(high.risk.density, bootstrap.hrisk.density, p.hrisk.target.cancer,
                "high risk kmers", element, strand, score, strand.label)
       # plot.me(low.risk.density, bootstrap.lrisk.density, p.lrisk.target.cancer,
       #         "low risk kmers", element, strand, score, strand.label)

      }, by = c(if (!both.strand) "strand", "element")]
    }
    dev.off()
  }

}
