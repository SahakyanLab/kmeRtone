#' Function performs an analysis of base composition including sequence 
#' frequency, PWM calculations, and G/C content at various window sizes.
#'
#' @param case A Coordinate class object or similar structure.
#' @param genome Genome class object or similar structure.
#' @param case.pattern String patterns to consider in the analysis.
#' @param output.path Output path for saving the analysis results.
#'
#' @importFrom data.table data.table rbindlist setnafill
#' @importFrom graphics barplot par plot lines legend
#' @importFrom grDevices cairo_pdf dev.off
#' @importFrom stringi stri_sub
#' @importFrom seqLogo makePWM seqLogo
#' @importFrom stats density
#' 
#' @export
countBaseComposition <- function(case, genome, case.pattern, output.path="./") {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (!is.null(output.path)) {
    dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
    cairo_pdf(paste0(output.path, "/exploration.pdf"), width = 10, height = 8,
              onefile = TRUE)
  }

  # options(scipen = 999)

  case$add_col_rep <- FALSE


  if (!is.null(case$single_len)) {

    combResults <- function(...) {

      seq.freq.dt <- rbindlist(lapply(list(...), `[[`, 1))[
        , .(N = sum(N)), by = seq]

      pwm.kmers <- Reduce("+", lapply(list(...), `[[`, 2))

      win.dt <- Reduce(function(dt1, dt2) dt1 + dt2, lapply(list(...),`[[`, 3))[
        , proportion := proportion / length(list(...))]

      return(list(seq.freq.dt, pwm.kmers, win.dt))
    }

    wins <- c(100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000,
              5000000)
    chr.lens <- genome$get_length(case$chr_names)
    wins <- wins[wins < (0.3 * min(chr.lens, na.rm = TRUE))]

    p <- progressr::progressor(along = coor$chr_names)

    results <- future.apply::future_lapply(coor$chr_names, function(chr.name) {

      p(chr.name)

      # ------------------------------------------------------------------------
      # Task 1: Get case sequence frequency

      case[chr.name]
      case$map_sequence(genome)
      seq.freq.dt <- case[chr.name][, .N, by = seq]

      # ------------------------------------------------------------------------
      # Task 2: Calculate PWM

      bases <- c("A", "C", "G", "T")
      pwm.win <- case$single_len + 100
      pwm.kmers <- matrix(0, nrow = 4, ncol = pwm.win, dimnames = list(bases))

      # Remove other than case.patterns
      case$coor[[chr.name]] <- case[chr.name][seq %in% case.pattern]

      # Expand to k-mer
      case[chr.name, state = "kmer", k = pwm.win]
      if (case$is_strand_sensitive)
        case[chr.name][strand == "-", start := - (start + pwm.win - 1)]

      # Calculate weight or probablity at every position
      for (i in 1:pwm.win) {
        base.tab <- case[chr.name][
          , {
            base <- stri_sub(genome[chr.name], abs(start + i - 1), length = 1)
            if (case$is_strand_sensitive) {
              base[strand == "-"] <- reverseComplement(base[strand == "-"])
            }
            table(base)
          }]
        base.tab <- base.tab[names(base.tab) %in% bases]
        pwm.kmers[names(base.tab), i] <- pwm.kmers[names(base.tab), i] +
          base.tab
      }

      # Return to original
      if (case$is_strand_sensitive)
        case[chr.name][strand == "-", start := - start - pwm.win + 1]
      case[chr.name, state = "case"]

      # ------------------------------------------------------------------------
      # Task 3: Calculate G and G|C content at various window size

      prop.interval <- 1000000
      win.dt <- data.table(proportion = 0:prop.interval / prop.interval,
                           key = "proportion")

      for (win in wins) {

        # Calculate G+C content in case
        win.GC <- countPointContext2(genome[chr.name][[1]],
                                     case[chr.name][, start],
                                     case$single_len, win, c("G", "C")) *
          if (case$is_strand_sensitive) 1 else 2
        win.GC.dt <- data.table(
          proportion = round(as.numeric(names(win.GC)) / win, digits = 6),
          N = win.GC,
          key = "proportion")
        win.dt[win.GC.dt, paste0("win", win, "_case_G|C") := N]

        # Calculate G content in case
        win.G <- c(
          countPointContext2(genome[chr.name][[1]],
                             case[chr.name][if (case$is_strand_sensitive)
                               strand == "+" else TRUE, start],
                             case$single_len, win, "G"),
          countPointContext2(genome[chr.name][[1]],
                             case[chr.name][if (case$is_strand_sensitive)
                               strand == "-" else TRUE, start],
                             case$single_len, win, "C"))
        win.G.dt <- data.table(
          proportion = round(as.numeric(names(win.G)) / win, digits = 6),
          N = win.G)[, .(N = sum(N)), by = proportion]
        setkey(win.G.dt, proportion)
        win.dt[win.G.dt, paste0("win", win, "_case_G") := N]

        # Calculate G+C content in genome
        win.GC <- countSlidingWindow2(genome[chr.name][[1]], win, c("G", "C"))
        win.GC.dt <- data.table(
          proportion = round(as.numeric(names(win.GC)) / win, digits = 6),
          N = win.GC * 2,
          key = "proportion"
        )
        win.dt[win.GC.dt, paste0("win", win, "_genome_G|C") := N]

        # Calculate G content in genome
        win.G <- c(countSlidingWindow2(genome[chr.name][[1]], win, "G"),
                   countSlidingWindow2(genome[chr.name][[1]], win, "C"))
        win.G.dt <- data.table(
          proportion = round(as.numeric(names(win.G)) / win, digits = 6),
          N = win.G
        )[, .(N = sum(N)), by = proportion]
        setkey(win.G.dt, proportion)
        win.dt[win.G.dt, paste0("win", win, "_genome_G") := N]

        # Calculate G and C content in all mid case pattern in the genome
        win.GC <- c(countMidPatternContext2(genome[chr.name][[1]], case.pattern,
                                            win, c("G", "C")),
                    countMidPatternContext2(genome[chr.name][[1]],
                                            reverseComplement(case.pattern),
                                            win, c("G", "C")))
        win.GC.dt <- data.table(
          proportion = round(as.numeric(names(win.GC)) / win, digits = 6),
          N = win.GC,
          key = "proportion"
        )
        win.dt[win.GC.dt, paste0("win", win, "_mid_G|C") := N]

        # Calculate G content in all mid case pattern in the genome
        win.G <- c(countMidPatternContext2(genome[chr.name][[1]], case.pattern,
                                           win, "G"),
                   countMidPatternContext2(genome[chr.name][[1]],
                                           reverseComplement(case.pattern),
                                           win, "C"))
        win.G.dt <- data.table(
          proportion = round(as.numeric(names(win.G)) / win, digits = 6),
          N = win.G,
          key = "proportion"
        )
        win.dt[win.G.dt, paste0("win", win, "_mid_G") := N]

        setnafill(win.dt, fill = 0)
      }

      return(list(seq.freq.dt, pwm.kmers, win.dt))
    }, future.seed = NULL)

    results <- combResults(results)

    par(oma = c(2, 2, 2, 2))

    # ------------------------------------------------------------------------
    # Task 1: Get case sequence frequency
    seq.freq.dt <- results[[1]]
    seq.freq.dt[
      ,
        barplot(N / sum(N) * 100, names.arg = seq,
                main = "Distribution of case pattern",
                xlab = "Case pattern",
                ylab = "Percentage")]

    # ------------------------------------------------------------------------
    # Task 2: Calculate PWM
    # Plot seqlogo at 100 window
    pwm.kmers <- results[[2]]
    pwm.kmers <- pwm.kmers / colSums(pwm.kmers)
    pwm.kmers <- seqLogo::makePWM(pwm.kmers)
    seqLogo::seqLogo(pwm.kmers)

    # ------------------------------------------------------------------------
    # Task 3: Calculate G and G|C content at various window size
    # Plot density of G and G|C content at various windows
    win.dt <- results[[3]]
    for (seq.content in c("G", "G|C")) {
      for (win in wins) {
        col.name <- paste0("win", win, "_case_", seq.content)
        dw.case <- win.dt[
          get(col.name) > 0,
          density(proportion, weights = get(col.name) / sum(get(col.name)))]
        col.name <- sub("case", "genome", col.name)
        dw.genome <- win.dt[
          get(col.name) > 0,
          density(proportion, weights = get(col.name) / sum(get(col.name)))]
        col.name <- sub("genome", "mid", col.name)
        dw.mid <- win.dt[
          get(col.name) > 0,
          density(proportion, weights = get(col.name) / sum(get(col.name)))]
        ymax <- max(dw.case$y, dw.genome$y, dw.mid$y)
        plot(dw.case, main = paste(seq.content, "content in", win, "window"),
             xlab = "Proportion", ylim = c(0, ymax), xlim = c(0, 1),
             col = "coral3", lwd = 3)
        lines(dw.genome, col = "cornflowerblue", lwd = 3)
        lines(dw.mid, col = "aquamarine4", lwd = 3)
        legend("topright", c("case", "genome", "possible case"),
               lty = 1, lwd = 3,
               col = c("coral3", "cornflowerblue", "aquamarine4"))
      }
    }
  }

  if (!is.null(output.path)) dev.off()
}
