#' Function performs an analysis of the distribution of genomic cases.
#' 
#' Check case distribution in replicates, chromosomes, and strands. 
#' Check case base composition and filter out other case.patterns. 
#' Then, it generates various plots like bar plots and Venn/Euler diagrams.
#'
#' @param case A Coordinate class object or similar structure for genomic data.
#' @param genome Genome class object or similar structure.
#' @param case.pattern String patterns to consider in the analysis.
#' @param output.path Output path for saving the analysis results.
#'
#' @importFrom data.table rbindlist setorder setorderv
#' @importFrom graphics barplot plot legend par mtext
#' @importFrom grDevices cairo_pdf dev.off
#' @importFrom stringi stri_sub
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor
#' @importFrom venneuler venneuler
#' 
#' @export
countDistribution <- function(case, genome, case.pattern, output.path="./") {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (!is.null(output.path)) {
    dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
    cairo_pdf(paste0(output.path, "/exploration.pdf"), width = 10, height = 8,
              onefile = TRUE)
  }

  case$add_col_rep <- TRUE

  p <- progressr::progressor(along = case$chr_names)

  # Counting
  count.dt <- future.apply::future_lapply(case$chr_names, function(chr.name) {

    p(chr.name)

    case[chr.name]
    if (!is.null(case.pattern)) case$map_sequence(genome)

    count.dt <- case[chr.name][
      if (!is.null(case.pattern)) seq %in% case.pattern else TRUE,
        .(chromosome = chr.name,
          N = if(!is.null(case$single_len)) .N * case$single_len
              else sum(end - start + 1)),
      by = eval(c("replicate", if (case$is_strand_sensitive) "strand"))]

    return(count.dt)
  }, future.seed = NULL) |> rbindlist()

  par(oma = c(2, 2, 2, 2), mar = c(7.1, 4.1, 4.1, 3.1))


  # Barplot
  count.dt[
    ,
    .(N = 100 * sum(N) / if (!is.null(case.pattern))
      genome$get_content(chromosome,
                         seq = if (case$is_strand_sensitive && strand == "-")
                           reverseComplement(case.pattern) else case.pattern)
      else genome$get_length(chromosome)),
    by = eval(c("chromosome", if (case$is_strand_sensitive) "strand"))][

    , barplot(matrix(N, nrow = 2),
              beside = TRUE,
              names.arg = unique(chromosome),
              ylab = "Percentage",
              col = c("lightgoldenrod4", "lightgoldenrodyellow"),
              main = "Case Distribution in Genome",
              las = 2)]
  mtext("Chromosome", side = 1, line = 4)

  if (case$is_strand_sensitive)
    legend(par("usr")[2], par("usr")[4], legend = count.dt[, unique(strand)],
           bty = "n", pch = 21, pt.cex = 1.2, cex = 0.8, horiz = FALSE,
           pt.bg = c("lightgoldenrod4", "lightgoldenrodyellow"),
           title = "strand", xpd = TRUE)

  # Venn or euler diagram
  if (length(unique(count.dt$replicate)) > 1) {
    setops <- count.dt[, .(N = sum(N)), by = replicate][
      , venneuler::venneuler(replicate, N)]

    cols <- c("cornflowerblue", "aquamarine3", "lightcoral", "mistyrose2",
              "khaki", "lightseagreen", "lemonchiffon1", "thistle2", "sienna3",
              "olivedrab3", "navajowhite2", "salmon2")
    alpha.val <- 0.5

    venn.cols <- cols[1:length(setops$labels)] |> addAlphaCol(alpha.val)
    names(venn.cols) <- setops$labels

    par(mar = c(2.1, 7.1, 4.1, 7.1), mgp = c(3, 1, 0))
    plot(setops, col.txt = NA, col = venn.cols, alpha = alpha.val,
         main = "Case Distribution in Replicates")

    setorderv(count.dt, "replicate")
    count.dt[grep("&", replicate), {
      reps <- strsplit(replicate, "&")[[1]]
      mixed.col <- mixColors(venn.cols[reps], alpha = alpha.val)
      names(mixed.col) <- replicate
      venn.cols <<- c(venn.cols, mixed.col)
      NULL
    }, by = replicate]

    perct.intsecs <- count.dt[, .(N = sum(N)), by = replicate][
      , N / sum(N) * 100] |> signif(3)
    names(perct.intsecs) <- count.dt$replicate |> unique()
    lab.intsecs <- paste0(perct.intsecs[names(venn.cols)], "% ",
                          gsub("&", " \U2229 ", names(venn.cols)))
    legend(par("usr")[1] - 0.1 * (par("usr")[2] - par("usr")[1]),
           par("usr")[4], legend = lab.intsecs,
           bty = "n", pch = 21, pt.bg = venn.cols, pt.cex = 1.5, cex = 0.8,
           horiz = FALSE, xpd = TRUE)
  }

  # Barplot for varied length
  if (is.null(case$single_len)) {

    len.dt <- lapply(case$chr_names, function(chr.name) {
      case[chr.name][, .(length = end - start + 1)][, .N, by = length]
    }) |> rbindlist()
    len.dt <- len.dt[, .(N = sum(N)), by = length]

    setorder(len.dt, length)
    len.dt[, barplot(N / sum(N) * 100, names.arg = length, xlab = "length",
                     ylab = "percentage", main = "Distribution of Case Length")]

  }

  if (!is.null(output.path)) dev.off()

}
