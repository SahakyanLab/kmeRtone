plotCancerGeneSusceptibilityDensity <- function(cancer.gene.census,
                                                transcript.susceptibility,
                                                tumour.exact.keyword,
                                                tumour.regex.keyword,
                                                output) {

  # Resolve regex keywords
  begin.rgx <- "((^|, )"
  end.rgx <- "($|,))"
  tumour.exact.keyword <- CJ(begin.rgx, tumour.exact.keyword, end.rgx)
  tumour.exact.keyword <- tumour.exact.keyword[, do.call(paste0, .SD)]
  tumour.keyword <- paste(c(tumour.exact.keyword, tumour.regex.keyword),
                          collapse = "|")

  # Change column name to match genome.annotation table for table merging
  setnames(cancer.gene.census, "Gene Symbol", "name2")

  # Add column "Tumour Types(Somatic)" and "Tumour Types(Germline)" to
  # transcript.susceptibility table
  transcript.susceptibility <-
    merge.data.table(transcript.susceptibility,
                     cancer.gene.census[, c("name2", "Tumour Types(Somatic)",
                                            "Tumour Types(Germline)")],
                     by = "name2", all.x = TRUE)


  # A helper function
  plot.me <- function(score, both.strand=FALSE) {

    # This command is specific to NCBIrefSeq table. I'm not sure about GENCODE
    setkey(transcript.susceptibility, strand, element, name, name2)

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

    if (both.strand) {

      idx.sense <- which(strand == "sense")
      idx.antisense <- which(strand == "antisense")

      base.content <- (base.content[idx.sense] +
                       base.content[idx.antisense]) / 2

      high.risk.content <- (high.risk.content[idx.sense] +
                            high.risk.content[idx.antisense]) / 2

      low.risk.content <- (low.risk.content[idx.sense] +
                           low.risk.content[idx.antisense]) / 2
    }

    # Total number of data. This is purely for both strand dataset because
    # we reduce the dataset by half in the if condition above. It won't affect
    # the single strand operation.
    n <- length(base.content)

    # Index of non-NaN number for high/low risk.
    # NaN as a result of 0/0 where there is no susceptible kmers at all.
    # Any mathematical operation involving NaN results in NaN.
    idx.hrisk.finite <- !is.nan(high.risk.content)[1:n]
    idx.lrisk.finite <- !is.nan(low.risk.content)[1:n]

    # Take index of all cancers and target cancer
    idx.all.cancer <- ((!is.na(`Tumour Types(Somatic)`)) |
                       (!is.na(`Tumour Types(Germline)`)))[1:n]
    idx.target.cancer <- (grepl(tumour.keyword, `Tumour Types(Somatic)`) |
                          grepl(tumour.keyword, `Tumour Types(Germline)`))[1:n]

    base.content.density <- list(
      density(base.content),
      density(base.content[idx.all.cancer]),
      density(base.content[idx.target.cancer])
    )

    high.risk.density <- list(
      density(high.risk.content[idx.hrisk.finite]),
      density(high.risk.content[idx.all.cancer & idx.hrisk.finite]),
      density(high.risk.content[idx.target.cancer & idx.hrisk.finite])
    )

    low.risk.density <- list(
      density(low.risk.content[idx.lrisk.finite]),
      density(low.risk.content[idx.all.cancer & idx.lrisk.finite]),
      density(low.risk.content[idx.target.cancer & idx.lrisk.finite])
    )

    # Plot parameters

    base.content.axes <- list(
      x = range(unlist(lapply(base.content.density, function(d) d$x))),
      y = range(unlist(lapply(base.content.density, function(d) d$y)))
    )

    high.risk.axes <- list(
      x = range(unlist(lapply(high.risk.density, function(d) d$x))),
      y = range(unlist(lapply(high.risk.density, function(d) d$y)))
    )

    low.risk.axes <- list(
      x = range(unlist(lapply(low.risk.density, function(d) d$x))),
      y = range(unlist(lapply(low.risk.density, function(d) d$y)))
    )

    legend.expr <- '
      legend("topright",
             legend = c("all genes", "all cancer genes",
                        "targeted cancer genes"),
             fill = c(c1, c2, c3)
      )'

    plot(base.content.density[[1]],
         main = paste("Distribution of base content in", element, "of",
                      c(if (length(strand) > 1) "both" else (strand)),
                      "strand"),
         xlim = base.content.axes$x,
         ylim = base.content.axes$y,
         col = c1 <- "black",
         xlab = "kmer fraction", lwd = 3
    )
    lines(base.content.density[[2]], col = c2 <- "orange", lwd = 3)
    lines(base.content.density[[3]], col = c3 <- "red", lwd = 3)
    eval(parse(text = legend.expr))

    plot(high.risk.density[[1]],
         main = paste("Distribution of high risk kmers in", element, "of",
                      c(if (length(strand) > 1) "both" else (strand)),
                      "strand"),
         xlim = high.risk.axes$x,
         ylim = high.risk.axes$y,
         col = c1,
         xlab = "kmer fraction", lwd = 3
    )
    lines(high.risk.density[[2]], col = c2, lwd = 3)
    lines(high.risk.density[[3]], col = c3, lwd = 3)
    eval(parse(text = legend.expr))

    plot(low.risk.density[[1]],
         main = paste("Distribution of low risk kmers in", element, "of",
                      c(if (length(strand) > 1) "both" else (strand)),
                      "strand"),
         xlim = low.risk.axes$x,
         ylim = low.risk.axes$y,
         col = c1,
         xlab = "kmer fraction", lwd = 3
    )
    lines(low.risk.density[[2]], col = c2, lwd = 3)
    lines(low.risk.density[[3]], col = c3, lwd = 3)
    eval(parse(text = legend.expr))

    NULL
    }, by = c(if (!both.strand) "strand", "element")]

  }

  scores <- names(transcript.susceptibility)[
                          grepl("risk|susceptible",
                                names(transcript.susceptibility))]
  scores <- unique(gsub("susceptible_kmers_|high_risk_kmers_|low_risk_kmers_",
                        "", scores))

  for (score in scores) {

    pdf(paste0(output, "/density_cancer_gene_susceptibility_", score, ".pdf"))

    plot.me(score)
    plot.me(score, both.strand = TRUE)

    dev.off()
  }

}
