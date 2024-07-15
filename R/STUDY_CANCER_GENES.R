#' Study k-mer composition of causal cancer genes from COSMIC Cancer Gene
#'    Census (CGC) database.
#'
#' Detail of Cancer Gene Census can be accessed and read at
#'    https://cancer.sanger.ac.uk/census
#'
#' @param cosmic.username COSMIC username i.e. registered email.
#' @param cosmic.password COSMIC password.
#' @param tumour.type.regex Regular expression for "Tumour Types" column in
#'    Cancer Gene Census table. Default is NULL.
#' @param tumour.type.exact Exact keywords for "Tumour Types" column in
#'    Cancer Gene Census table. Default is NULL.
#' @param cell.type Type of cell: "somatic" or "germline". Default is "somatic".
#' @param genic.elements.counts.dt Genic element count table generated from
#'    STUDY_GENIC_ELEMENTS.
#' @param output.dir A directory for the outputs.
#'
#' @return An output directory containing plots.
#'
#' @importFrom data.table fread fwrite setorder setnames rbindlist
#' @importFrom stringi stri_sub stri_count_regex stri_length stri_paste stri_sub_replace_all stri_replace_all_regex stri_split_fixed
#' @importFrom Biostrings reverseComplement
#' @importFrom graphics layout boxplot plot points legend arrows axis mtext barplot abline box grconvertX grconvertY
#'    matlines matplot polygon rasterImage rect text title
#' @importFrom grDevices pdf png
#' @importFrom stats C end median runif sd start wilcox.test
#' 
#' @export
STUDY_CANCER_GENES <- function(cosmic.username, cosmic.password,
                               tumour.type.regex=NULL, tumour.type.exact=NULL,
                               cell.type="somatic", genic.elements.counts.dt,
                               output.dir="study_cancer_genes/") {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # Convert exact keywords to regex
  rgx.exact <- "^WORD$|^WORD,|, WORD,|, WORD$"
  tumour.type.regex <- paste(
    c(tumour.type.regex,
      sapply(tumour.type.exact, function(keyword) gsub("WORD", keyword,
                                                     rgx.exact))),
  collapse = "|")

  if (is.character(genic.elements.counts.dt)) {
    counts <- fread(genic.elements.counts.dt, showProgress = FALSE)
  } else {
    counts <- genic.elements.counts.dt
  }

  # Get cancer genes -----------------------------------------------------------
  cgc <- getCOSMICcancerGeneCensus(email = cosmic.username,
                                   password = cosmic.password)

  # Resolve column to filter i.e. either somatic or germline
  filter.col <- paste0("Tumour Types(", stri_trans_totitle(cell.type), ")")

  # Filter cancer genes based on regex
  sel.idx <- cgc[, get(filter.col) %like% tumour.type.regex]

  sel.cancer.genes <- cgc[sel.idx, `Gene Symbol`]
  oth.cancer.genes <- cgc[!sel.idx, `Gene Symbol`]

  cgc[`Gene Symbol` %in% sel.cancer.genes, selected := TRUE]
  cgc[is.na(selected), selected := FALSE]

  fwrite(cgc, paste0(output.dir, "/cancer_gene_census.csv"),
         showProgress = FALSE)

  # Mark cancer genes in the count table
  gene.name.col <- names(counts)[names(counts) %in% c("name2", "geneName")]
  counts[get(gene.name.col) %in% sel.cancer.genes, cancer := "selected"]
  counts[get(gene.name.col) %in% oth.cancer.genes, cancer := "other"]

  fwrite(counts, paste0(output.dir, "/genic_elements_counts_cancer.csv"),
         showProgress = FALSE)

  # Check if any cancer genes in CGS is not in genePred
  missing.genes <- counts[
    cancer %in% c("selected", "other"),
    sel.cancer.genes[!sel.cancer.genes %in% unique(get(gene.name.col))]
  ]
  if (length(missing.genes) > 0)
    warning("The following cancer genes are not found in genePred: ",
            missing.genes)

  # Calculate density ----------------------------------------------------------
  # Density of selected/other cancer genes vs. all other none cancer genes.
  # Since number of all other none cancer genes is too high compare to cancer
  # genes, repeated sub-sampling of all other none cancer genes is performed to
  # match the number of the selected/other cancer genes.

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
    top_kmer = sum(top_kmer),
    cancer = cancer[1]
  ), by = c(gene.name.col, "element")][, `:=`(
    "susceptible k-mer content" = susceptible_kmer / total_kmer,
    "highly susceptible k-mer content" = top_kmer / susceptible_kmer
  )]
  # NaN because zero k-mer content. 0 / 0 = NaN
  counts.both[is.nan(`highly susceptible k-mer content`),
              "highly susceptible k-mer content" := 0]

  # Order table by element appearance in the plot.
  element.sort <- 1:6
  names(element.sort) <- c("upstream", "CDS", "3'-UTR", "5'-UTR", "intron",
                           "downstream")
  counts[, rank := element.sort[element]]
  setorder(counts, rank, -strand)
  counts[, rank := NULL]
  counts.both[, rank := element.sort[element]]
  setorder(counts.both, rank)
  counts.both[, rank := NULL]

  for (cancer.grp in c("selected", "other")) {

    n.cancer.genes <- ifelse(cancer.grp == "selected", length(sel.cancer.genes),
                             length(oth.cancer.genes))

    pdf(paste0(output.dir, "/", cancer.grp, "_cancer_genes_density_plots.pdf"),
        paper = "a4r", width = 11.69, height = 8.27)

    # Plot layout.
    # Panel 1,3,5,7,9,11 are density plots.
    # Panel 2,4,6,8,10,12 are significance score plots.
    # Panel 13 is title.
    # Panel 14 is gene construct drawing.
    # Panel 15 is ylab for density plots.
    # Panel 16 is ylab for significance score plots.
    # Panel 17 is xlab for all plots.
    lay.mat <- matrix(c(0, rep(13, 3), 0, seq(2, 6, 2), 0, seq(1, 6, 2), 0,
                        rep(14, 3), 16, seq(8, 12, 2), 15, seq(7, 12, 2), 0, 17,
                        0, 0),
                      nrow = 7, byrow = TRUE)
    layout(lay.mat, heights = c(2, 3, 5, 3, 3, 5, 0.875),
          widths = c(0.5, 4, 4, 4))

    boot.num <- 1e5
    subsample.non.cancer.genes <- counts[element != "IGR" & strand == "sense" &
      is.na(cancer),
      replicate(boot.num, sample(get(gene.name.col), n.cancer.genes),
                simplify = TRUE)]

    for (kmer.content in c("susceptible k-mer content",
                           "highly susceptible k-mer content")) {

      for (strand.type in c("both", "sense", "antisense")) {

        (if (strand.type == "both") counts.both else counts)[element != "IGR" &
                if (strand.type == "both") TRUE else strand == strand.type, {

          # Initially calculate density of all genes to get k-mer proportion
          # reference. Outliers are not removed because we could get
          # distribution around outliers for cancer genes.
          kmer.content.ref <- density(get(kmer.content), na.rm = TRUE)$x

          density.ref <- function(kmer.content) {
            density(kmer.content,
                    n = length(kmer.content.ref),
                    from = kmer.content.ref[1],
                    to = kmer.content.ref[length(kmer.content.ref)],
                    na.rm = TRUE)$y
          }

          p <- progressr::progressor(steps = boot.num)
          t1 <- Sys.time()

          # Bootstrap genes and calculate density - boot.num (100k) times
          dens.none.cancer.genes.mat <- apply(subsample.non.cancer.genes, 2,
            function(gene.names) {
              p(paste("Bootstraping:", kmer.content, "in", element, "of",
                      strand.type, "strands"))
              density.ref(get(kmer.content)[get(gene.name.col) %in% gene.names])
            })

          time.taken <- Sys.time() - t1
          message(paste("Bootstraping genes for ", kmer.content, " in ", element, " of ",
              strand.type, " strands in ", cancer.grp, " cancer genes...",
              round(time.taken, 2), " ", attr(time.taken, "units"), "\n",
              sep = ""))

          dens.cancer.genes <- density.ref(get(kmer.content)[cancer ==
            cancer.grp])

          dens.none.cancer.genes <- density.ref(get(kmer.content)[
            is.na(cancer)])

          # Perform wilcox.test
          sig.scores <- sapply(seq_along(kmer.content.ref), function(i) {

            # Calculate p-value
            side <- ifelse(dens.cancer.genes[i] >
                            median(dens.none.cancer.genes.mat[i, ]),
                          "greater", "less")
            p.val <- wilcox.test(x = dens.cancer.genes[i],
                                y = dens.none.cancer.genes.mat[i, ],
                                alternative = side)$p.value

            # Calculate significance score
            sig.score <- -log(p.val)
            if (side == "less") sig.score <- sig.score * -1

            return(sig.score)
          })

          # Plotting -----------------------------------------------------------
          par(mgp = c(3, 0.7, 0)) # Move tick labels close to the tick.

          # Panel 1,3,5,7,9,11: density plots
          par(mar = c(1, 1, 0, 1))
          par.mar <- par("mar")
          dens.range <- range(dens.none.cancer.genes.mat, dens.cancer.genes,
                              dens.none.cancer.genes)

          # Initialise empty plot.
          plot(NULL, xlim = range(kmer.content.ref), ylim = dens.range,
              axes = FALSE, ann = FALSE)
          axis(1, col = NA, col.ticks = "azure3")
          axis(2, col = NA, col.ticks = "azure3")
          legend("topright", element, bty = "n", cex = 0.8,
                 text.col = "azure3")

          # Rasterize first to avoid massive amount of line vectors
          # Get X and Y coordinates in inches
          par.usr <- par("usr")
          g.x <- grconvertX(par.usr[1:2], "user", "inches")
          g.y <- grconvertY(par.usr[3:4], "user", "inches")
          g.width <- max(g.x) - min(g.x)
          g.height <- max(g.y) - min(g.y)

          tmp.png <- tempfile(fileext = ".png")
          png(tmp.png, width = g.width, height = g.height, units = "in",
              res = 300, bg = "transparent")
          par(mar = rep(0, 4))
          plot(NULL, xlim = par.usr[1:2], ylim = par.usr[3:4],
              xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)

          matlines(x = kmer.content.ref, y = dens.none.cancer.genes.mat,
                  col = addAlphaCol("gray", 1 / 200), type = "l", lty = 1)

                    dev.off() # Rasterize ends here.

          par(mar = par.mar)

          # Read the rasterized lines and shade.
          rasterImage(image = png::readPNG(tmp.png),
                      xleft = par.usr[1],
                      ybottom = par.usr[3],
                      xright = par.usr[2],
                      ytop = par.usr[4])

          # Shade significantly different density region
          sig.up <- R.utils::seqToIntervals(which(sig.scores > -log(0.05)))
          sig.down <- R.utils::seqToIntervals(which(sig.scores < log(0.05)))

          if (length(sig.up) > 0) {
            rect(xleft = kmer.content.ref[sig.up[, "from"]],
                xright = kmer.content.ref[sig.up[, "to"]],
                ybottom = par.usr[3], ytop = par.usr[4],
                col = addAlphaCol("indianred1", 0.2),
                border = NA)
          }

          if (length(sig.down) > 0) {
            rect(xleft = kmer.content.ref[sig.down[, "from"]],
                xright = kmer.content.ref[sig.down[, "to"]],
                ybottom = par.usr[3], ytop = par.usr[4],
                col = addAlphaCol("darkseagreen3", 0.2),
                border = NA)
          }

          lines(x = kmer.content.ref, y = dens.none.cancer.genes,
                col = "black")

          lines(x = kmer.content.ref, y = dens.cancer.genes,
                col = "orangered4")

          box(col = "azure3")

          # Panel 2,4,6,8,10,12 - Significance plots
          # Plot significance score above the density plot
          par(mar = c(1, 1, 1, 1))
          plot(NULL, xlim = range(kmer.content.ref),
              ylim = c(min(-4, sig.scores), max(4, sig.scores)),
              axes = FALSE, ann = FALSE, xaxs = 'i')
          axis(2, col = NA, col.ticks = "azure3")

          # Shade
          rect(xleft = par("usr")[1],
              xright = par("usr")[2],
              ybottom = -log(0.05),
              ytop = par("usr")[4],
              col = addAlphaCol("indianred1", 0.2),
              border = NA)

          rect(xleft = par("usr")[1],
              xright = par("usr")[2],
              ybottom = par("usr")[3],
              ytop = log(0.05),
              col = addAlphaCol("darkseagreen3", 0.2),
              border = NA)

          abline(h = 0, lty = 2, col = "grey")

          lines(x = kmer.content.ref, y = sig.scores)

          box(col = "azure3")

        }, by = element]

        # Panel 13 - Title
        par(mar = c(0, 1, 0, 1))
        plot(1, type = "n", axes = FALSE, ann = FALSE)
        text(x = par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2,
            y = par("usr")[3] + (par("usr")[4] - par("usr")[3]) / 2,
            labels = paste0(kmer.content, " in ", strand.type, " strand"),
            cex = 1.5, font = 4)

        # Panel 14 - Gene construct drawing
        par(mar = c(0, 1, 2, 1))
        plot(NULL, xlim = c(0, 11), ylim = c(0, 10), xaxs = 'i', yaxs = 'i',
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

        # Panel 15 - ylab for density plots.
        par(mar = c(1, 1, 0, 0))
        plot(1, type = "n", axes = FALSE, ann = FALSE)
        mtext("density", side = 4, line = -2.5, cex = 0.8)

        # Panel 16 - ylab for significance score plots.
        par(mar = c(1, 1, 1, 0))
        plot(1, type = "n", axes = FALSE, ann = FALSE)
        mtext("sig. score", side = 4, line = -2.5, cex = 0.8)

        # Panel 17 - xlab for all plots.
        par(mar = c(1, 1, 0, 1))
        plot(1, type = "n", axes = FALSE, ann = FALSE)
        mtext("k-mer content", side = 3, line = -2.5, cex = 0.8)

      }
    }
    dev.off()
    NULL
  }

}