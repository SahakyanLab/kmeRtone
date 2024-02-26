#' Study k-mer composition across species.
#'
#' Analysis of distribution of highly enriched k-mers across species.
#'
#' @param kmer.table A data.table of kmer table or path to it.
#' @param kmer.cutoff Percentage of extreme kmers to study. Default to 5
#'    percent.
#' @param selected.extremophiles A vector of selected extremophile species. e.g.
#'    c("Deinococcus soli", "Deinococcus deserti")
#'    The best representative will be selected from the assembly summary.
#' @param other.extremophiles A vector of other extremophile species. These are
#'    used as a control to compare with the selected extremophiles.
#' @param k K-mer size.
#' @param central.pattern K-mer's central patterns. Default is NULL.
#' @param output.dir A directory for the outputs.
#'
#' @return An output directory containing plots.
#'
#' @export
STUDY_ACROSS_SPECIES <- function(kmer.table, kmer.cutoff=5, k,
                                 central.pattern=NULL, selected.extremophiles,
                                 other.extremophiles,
                                 output.dir="study_across_species/") {

  if (is.character(kmer.table))
    kmer.table <- fread(kmer.table, showProgress = FALSE)

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  # Select dsDNA genomes -------------------------------------------------------
  organism.groups <- c("archaea", "bacteria", "fungi", "invertebrate", "plant",
                       "protozoa", "vertebrate_mammalian", "vertebrate_other",
                       "viral")
  asm <- rbindlist(lapply(organism.groups, function(organism.group) {
    selectGenomesForCrossSpeciesStudy(organism.group = organism.group,
                                      db = "refseq")
  }))

  # Select only dsDNA virus
  vmr <- getICTVvirusMetadataResource()[Genome.composition == "dsDNA"]
  asm[organism_group == "viral",
      idx := organism_name %in% vmr[, c(`Virus.name(s)`, Species)]]
  asm[organism_group != "viral", idx := TRUE]
  asm <- asm[(idx)][, idx := NULL]

  # Select extremophile genomes ------------------------------------------------
  asm.refseq <- getNCBIassemblySummary(organism.group = "all", db = "refseq")
  asm.genbank <- getNCBIassemblySummary(organism.group = "all", db = "genbank")
  # selected.extremophiles <- c(
  #     "Deinococcus radiodurans", "Rubrobacter xylanophilus",
  #     "Deinococcus deserti", "Deinococcus gobiensis",
  #     "Rubrobacter radiotolerans", "Deinococcus soli",
  #     "Cellulosimicrobium cellulans", "Deinococcus wulumuqiensis")
  # other.extremophiles = c(
  #     "Deinococcus maricopensis", "Deinococcus proteolyticus",
  #     "Deinococcus geothermalis", "Chroococcidiopsis thermalis",
  #     "Deinococcus peraridilitoris", "Hymenobacter glacialis",
  #     "Colwellia marinimaniae", "Moritella marina", "Streptomyces radiopugnans",
  #     "Hymenobacter roseosalivarius")
  tardigrade = c("Hypsibius exemplaris", "Ramazzottius varieornatus")
  extremophile.names <- list(selected_extremophile = selected.extremophiles,
                             other_extremophile = other.extremophiles,
                             tardigrade = tardigrade)
  extreme.asm <- rbindlist(
    lapply(seq_along(extremophile.names), function(i) {
      rgx <- paste0("(", paste(extremophile.names[[i]], collapse = ")|("), ")")
      asm <- asm.refseq[organism_name %like% rgx] |>
        selectRepresentativeFromASM()
      asm[, organism_group := names(extremophile.names)[i]]

      # Check if any extremophile not found, use genbank database
      if (asm[, .N] < length(extremophile.names[[i]])) {
        miss.extr <- sapply(extremophile.names[[i]], function(extr) {
          asm[organism_name %like% extr, .N] == 0
        })
        miss.extr <- extremophile.names[[i]][miss.extr]
        rgx <- paste0("(", paste(miss.extr, collapse = ")|("), ")")
        asm <- rbind(asm,
                     selectRepresentativeFromASM(
                       asm.genbank[organism_name %like% rgx])[
                         , organism_group := names(extremophile.names)[i]])
      }

      # Check again. If still missing, report to user.
      if (asm[, .N] < length(extremophile.names[[i]])) {
        miss.extr <- sapply(extremophile.names[[i]], function(extr) {
          asm[organism_name %like% extr, .N] == 0
        })
        miss.extr <- extremophile.names[[i]][miss.extr]
        cat(paste(miss.extr[-length(miss.extr)], collapse = ", "), "and",
            miss.extr[length(miss.extr)], " are not found in both  NCBI refseq",
            " and genbank database.\n")
        response <- readline("Continue without these extremophiles? (y/n): ")
        if (response == "n") stop("Stop R")
      }

      return(asm)
    })
  )

  asm <- rbind(asm, extreme.asm)

  # Save assembly summary for future reference.
  fwrite(asm, paste0(output.dir, "/assembly_summary.csv"),
         showProgress = FALSE)

  # Remove report_summary, if any
  unlink(paste0(output.dir, "/assembly_reports.csv"))

  # Count k-mer composition ----------------------------------------------------
  setorder(kmer.table, -z)
  extreme.kmers <- kmer.table[1:round(kmer.cutoff / 100 * .N), kmer]
  small.genomes <- c("bacteria", "viral", "archaea", "protozoa", "fungi")
  central.pattern.size <- unique(nchar(central.pattern))
  if (length(central.pattern.size) > 1)
    stop("Central pattern size should be the same.")

  p <- progressr::progressor(steps = nrow(asm))

  asm[, `:=`(top_kmer_count = {

    #cat(assembly_accession, organism_name, organism_group, "\n")
    p(paste(assembly_accession, organism_name, organism_group))

    genome <- loadGenome(genome.name = assembly_accession,
                         fasta.style = "NCBI",
                         mask = "none")

    report <- genome$get_assembly_report()

    # Filter to only include chromosome sequence
    # Assembly report is ordered by chromosome/plasmid;
    # the chromosomes/plasmids are followed by unlocalized scaffolds.
    # Unplaced scaffolds are listed at the end.
    # So just take the first row as reference.
    report <- report[`Sequence-Role` == `Sequence-Role`[1] &
                       `Assigned-Molecule-Location/Type` ==
                       `Assigned-Molecule-Location/Type`[1]]

    # Save filtered assembly report as a reference.
    fwrite(report, paste0(output.dir, "/assembly_reports.csv"),
           append = TRUE, showProgress = FALSE)

    sel.chrs <- genome$avail_seqs[genome$avail_seqs %in%
                                    report[, c(`RefSeq-Accn`, `GenBank-Accn`)]]

    total.kmer.cnt <- 0
    susceptible.kmer.cnt <- 0
    top.kmer.cnt <- 0

    for (sel.chr in sel.chrs) {

      # kmers <- buildKmerTable(genome[sel.chr], k, method = "auto",
      #                         remove.N = FALSE)
      kmers <- buildMidPatternKmerTable(
        dna.seqs = genome[sel.chr], k,
        mid.patterns = c(central.pattern, reverseComplement(central.pattern)),
        remove.N = FALSE)

      top.kmer.cnt <-
        sum(top.kmer.cnt,
            kmers[kmer %in% c(extreme.kmers,
                              reverseComplement(extreme.kmers)), sum(N)])

      # susceptible.kmer.cnt <-
      #   sum(susceptible.kmer.cnt,
      #       kmers[stri_sub(kmer, from = (k - central.pattern.size) / 2 + 1,
      #                      length = central.pattern.size) %in%
      #               c(central.pattern, reverseComplement(central.pattern)),
      #             sum(N)])
      susceptible.kmer.cnt <- sum(susceptible.kmer.cnt, kmers[, sum(N)])

      total.kmer.cnt <- sum(total.kmer.cnt, sum(genome$seq_len - k + 1) * 2)

    }
    if (total.kmer.cnt == 0) {
      stop(assembly_accession, " has zero total k-mer count.")
    }
    top.kmer.cnt
  },
  susceptible_kmer_count = susceptible.kmer.cnt,
  total_kmer_count = total.kmer.cnt),
  by = seq_len(nrow(asm))]

  # Save assembly summary table with count.
  fwrite(asm, paste0(output.dir, "/assembly_summary.csv"),
         showProgress = FALSE)

  # Plot k-mer composition -----------------------------------------------------
  # Highly susceptible k-mer content vs. Susceptible k-mer content
  pdf(file = paste0(output.dir, "/plots.pdf"), paper = "a4r",
      width = 11.69, height = 8.27)

  org.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#FFFF33", "#A65628", "#F781BF", "#999999") |> addAlphaCol(0.5)
  names(org.cols) <- organism.groups

  asm[, `:=`("susceptible k-mer content" = susceptible_kmer_count /
               total_kmer_count,
             "highly susceptible k-mer content" = top_kmer_count /
               susceptible_kmer_count)]

  # Plot layout
  # Panel 1 - scatter plot
  # Panel 2 - density plot for susceptible k-mer content
  # Panel 3 - density plot for highly susceptible k-mer content
  # Panel 4 - legend
  lay.mat <- matrix(c(2, 4,
                      1, 3),
                    nrow = 2, byrow = TRUE)
  layout(lay.mat, heights = c(2, 5), widths = c(5, 2))

  # Panel 1 - scatter plot
  par(mar = c(5.1, 4.1, 0, 0))

  # Initialise an empty plot
  plot(`highly susceptible k-mer content` ~ `susceptible k-mer content`,
       data = asm, type = "n",
       xlab = "susceptible k-mer content",
       ylab = "highyly susceptible k-mer content")

  # Rasterize first to avoid massive amount of pch vectors
  # Get X and Y coordinates in inches
  par.usr <- par("usr")
  par.mar <- par("mar")
  g.x <- grconvertX(par.usr[1:2], "user", "inches")
  g.y <- grconvertY(par.usr[3:4], "user", "inches")
  g.width <- max(g.x) - min(g.x)
  g.height <- max(g.y) - min(g.y)

  tmp.png <- tempfile(fileext = ".png")
  png(tmp.png, width = g.width, height = g.height, units = "in",
      res = 300, bg = "transparent")
  par(mar = rep(0, 4))
  plot(NULL, xlim = par.usr[1:2], ylim = par.usr[3:4],
       xaxs = 'i', yaxs = 'i', axes = FALSE, ann = FALSE)

  # Order by organism_group number
  org.order <- seq_along(asm[, unique(organism_group)])
  names(org.order) <- asm[, .N, by = organism_group][order(-N)]$organism_group
  asm[, temp := org.order[organism_group]]
  setorder(asm, temp)[, temp := NULL]

  for (organism.group in asm[, unique(organism_group)])
    asm[organism_group == organism.group, {

      points(x = `susceptible k-mer content`,
             y = `highly susceptible k-mer content`,
             pch = 16,
             col = org.cols[organism_group])

      NULL
    }]

  dev.off() # Rasterisation ends here.

  par(mar = par.mar)

  # Read the rasterized lines and shade.
  rasterImage(image = png::readPNG(tmp.png),
              xleft = par.usr[1],
              ybottom = par.usr[3],
              xright = par.usr[2],
              ytop = par.usr[4])

  greek.letters <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta",
                     "eta", "theta", "iota", "kappa", "lambda", "mu", "nu",
                     "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon",
                     "phi", "chi", "psi", "omega")
  asm[organism_group == "tardigrade", label := greek.letters[1:.N]]
  asm[organism_group == "selected_extremophile", label := LETTERS[1:.N]]
  asm[organism_group == "other_extremophile", label := letters[1:.N]]

  # Save assembly summary table with labels on plot.
  fwrite(asm, paste0(output.dir, "/assembly_summary.csv"),
         showProgress = FALSE)

  # Label extremophiles
  asm[organism_group %in% c("tardigrade", "selected_extremophile",
                            "other_extremophile"),
      {

        text(x = `susceptible k-mer content`,
             y = `highly susceptible k-mer content`,
             labels = sapply(label, function(l)
               eval(parse(text = paste0("expression(", l, ")")))),
             cex = 0.8, family = "mono")

        NULL
      }
  ]

  # Panel 2 - density plot for susceptible k-mer content
  #par(mar = c(0.5, 4.1, 1, 0))
  par(mar = c(0, 4.1, 1, 0))
  dens <- lapply(organism.groups, function(g) {
    density(asm[organism_group == g][["susceptible k-mer content"]])
  })
  Xs <- sapply(dens, `[[`, "x")
  Ys <- sapply(dens, `[[`, "y")

  matplot(Xs, Ys, type = "l", lty = 1, col = org.cols, lwd = 2, ann = FALSE,
          axes = FALSE, xlim = par.usr[1:2], xaxs = 'i')

  # Panel 3 - density plot for highly susceptible k-mer content
  #par(mar = c(5.1, 0.5, 0, 1))
  par(mar = c(5.1, 0, 0, 1))
  dens <- lapply(organism.groups, function(g) {
    density(asm[organism_group == g][["highly susceptible k-mer content"]])
  })
  Xs <- sapply(dens, `[[`, "x")
  Ys <- sapply(dens, `[[`, "y")

  matplot(Ys, Xs, type = "l", lty = 1, col = org.cols, lwd = 2, ann = FALSE,
          axes = FALSE, ylim = par.usr[3:4], yaxs = 'i')

  # Panel 4 - legend
  # Legend 1 - Prokaryote: bacteria, archaea
  # Legend 2 - Eukaryote: protozoa, fungi, vertebrate_mammalian,
  #                       vertebrate_other, invertebrate, plant
  # Legend 3 - Virus
  # Legend 4 - Extremophiles
  prokaryotes <- c("bacteria", "archaea")
  eukaryotes <- c("protozoa", "fungi", "plant", "vertebrate_mammalian",
                  "vertebrate_other", "invertebrate")
  par(mar = c(0, 0, 1, 1))
  plot(1, type = "n", axes = FALSE, ann = FALSE)
  space <- 0.05 * (par("usr")[4] - par("usr")[3])
  lgnd1 <- legend(x = par("usr")[1],
                  y = par("usr")[4],
                  legend = gsub("_", " ", names(org.cols[prokaryotes])),
                  col = org.cols[prokaryotes],
                  pch = 16, bty = 'n', ncol = 2, text.width = NA,
                  title = "  Prokaryote", title.adj = 0,
                  y.intersp = 0.8)
  lgnd2 <- legend(x = par("usr")[1],
                  y = min(lgnd1$text$y) - space,
                  legend = gsub("_", " ", names(org.cols[eukaryotes])),
                  col = org.cols[eukaryotes],
                  pch = 16, bty = 'n', ncol = 2, text.width = NA,
                  title = "  Eukaryote", title.adj = 0,
                  y.intersp = 0.8)
  lgnd3 <- legend(x = par("usr")[1],
                  y = min(lgnd2$text$y) - space,
                  legend = "virus",
                  col = org.cols["viral"],
                  pch = 16, bty = 'n', ncol = 2, text.width = NA,
                  y.intersp = 0.8)
  etxm.labs <- asm[organism_group %in% c("tardigrade", "selected_extremophile",
                                         "other_extremophile"),
                   paste0("expression(", label[1], "*`-`*", label[.N], ")"),
                   by = organism_group][match(organism_group,
                                              c("selected_extremophile",
                                                "other_extremophile",
                                                "tardigrade"))]
  par(family = "mono")
  lgnd4a <- legend(x = par("usr")[1],
                   y = min(lgnd3$text$y) - space,
                   legend = sapply(etxm.labs$V1,
                                   function(ex) eval(parse(text = ex))),
                   x.intersp = 0, cex = 0.85, y.intersp = 1.17, bty = "n")
  par(family = "")
  lgnd4b <- legend(x = lgnd4a$rect$left + lgnd4a$rect$w,
                   y = min(lgnd3$text$y) - space,
                   legend = gsub("_", " ", etxm.labs$organism_group),
                   x.intersp = -0.9, bty = "n", yjust = 0.95)

  rect(xleft = lgnd1$rect$left,
       ybottom = lgnd4b$rect$top - lgnd4b$rect$h,
       xright = max(lgnd1$rect$left + lgnd1$rect$w,
                    lgnd2$rect$left + lgnd2$rect$w,
                    lgnd3$rect$left + lgnd3$rect$w,
                    lgnd4b$rect$left + lgnd4b$rect$w),
       ytop = lgnd1$rect$top)

  dev.off()
}
