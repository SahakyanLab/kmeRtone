drawGene <- function(genomic.coordinate, colours=NULL) {

  elements <- c("upstream", "UTR5", "CDS", "intron", "UTR3", "downstream",
                "UTR")

  # Set element colours
  if (is.null(colours)) {
    set.seed(3)
    colours <- sample(colours(), length(elements))
  }
  names(colours) <- elements

  dt <- genomic.coordinate[element %in% elements]

  gene.name <- dt[, unique(name)]

  # Update elements. Some elemets can be not in existence
  exist.elements <- dt[, unique(element)]
  elements <- elements[elements %in% exist.elements]
  colours <- colours[elements]

  # locate any overlaps among damage and control regions
  setkey(dt, start, end)
  dt <- dt[, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < 0))]

  #Assign bin number for overlaps
  assignBin <- function(start, end) {

    N <- length(start)

    if (N == 1) {
      bin <- 1
    } else {

      bin <- c(1, rep(NA, N-1))

      for (i in 2:N) {

        end.before <- end[1:(i-1)]
        bin.before <- bin[1:(i-1)]

        if (sum(start[i] <= end.before) >= max(bin.before)) {

          bin[i] <- max(bin.before) + 1

        } else if (sum(start[i] <= end.before) == 0) {

          bin[i] <- 1

        } else {

          # get bin number for overlapping region
          bin.overlap <- bin.before[start[i] <= end.before]

          # assign bin number that does not belong to the overlaps
          bin[i] <- bin.before[!bin.before %in% bin.overlap][1]

        }
      }
    }
    return(bin)
  }

  # assign bin number for every continuous regions i.e. non-overlapping
  dt[, bin := assignBin(start, end), by = .(name, group)]

  height <- 0.2
  xlim <- c(min(dt$start), max(dt$end))
  bins <- dt$bin # to detect overlap and draw overlap above
  plot.new()
  plot.window(xlim = xlim, c(0, max(bins)*length(gene.name) * (height + 0.5)))
  ybottom <- bins * (0.2 + height) - height

  ybottoms <- lapply(elements, function(element_) {
                ybottom[dt$element == element_]})
  names(ybottoms) <- elements

  i=0
  for (name_ in dt[, unique(name)]){

    ybottoms_ <- NULL
    for (element_ in elements){

      start <- dt[element == element_ & name == name_, start]
      end <- dt[element == element_ & name == name_, end]
      ybottom <- ybottoms[[element_]][dt[element == element_, name == name_]]
      ybottoms_ <- c(ybottoms_, ybottom)

      rect(start - 0.5, ybottom+i, end + 0.5,
           ybottom+i + height, col = colours[[element_]], border=NA)

    }
    text(x=dt[name == name_, min(start)] - 0.5, y=min(ybottoms_)+i-0.05,
         labels = name_, adj = 0, cex = 0.8)

    i <- i+0.5
  }

  par(xpd=TRUE)

  legend("topright", legend=elements,
         fill=colours, cex=0.8, bty = "n")

  #title(xlab=paste("Transcript", paste( dt[, unique(name)], collapse = ", ")))

  axis(1)

}
