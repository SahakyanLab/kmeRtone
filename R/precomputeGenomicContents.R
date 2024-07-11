#' Function calculates various genomic content metrics based on the provided genome object.
#'
#' @param genome An object of class 'NCBI_Genome' containing genomic information.
#'
#' @return A data.table containing calculated genomic content metrics.
#'
#' @importFrom data.table fread fwrite set
#'
#' @export
loadGenomicContents <- function(genome) {

  # Support
  # Contents: "G" and "G|C"
  # Windows: 100, 500, 1000, 5000, 10000, 50000, 100000, 200000
  # Middle Seq: G, C, A, T, CT, CC, TC, TT
  # Columns: proportion, win100k_G, win1000k_G|C, win100000k_G_CT, etc.
  # No mitochondria genome

  genome$get_length("")
  chr.names <- fread(genome$info_file)$chromosome |> unique()
  chr.names <- chr.names[!grepl("chrM", chr.names)]

  win.file <- paste0(genome$root_path, "/window.csv")
  if (file.exists(win.file)) {
    win.dt <- fread(win.file)
    return(win.dt)
  } else {
    win.dt <- data.table(proportion = round(0:1000000 / 1000000, digits = 6))
  }

  wins <- c(100, 500, 1000, 5000, 10000, 50000, 100000, 200000)
  for (win in wins) {

    message(paste("At window", win, "\n"))

    # Calculate G content in genome on window of both strand
    win.dt[, paste0("win", win, "_G") :=
      genome$count_win_content(chr.names, c("G", "C"), win)[
      , .(proportion, G = G + C)]$G]

    mid.seqs <- c("G", "A", "CC", "CT", "TC", "TT")
    mid.seqs <- c(mid.seqs, reverseComplement(mid.seqs))
    for (mid.seq in mid.seqs) {
      set(win.dt, j = paste0("win", win, "_G_", mid.seq),
          value = genome$count_win_content(chr.names, c("G", "C"), win,
                                                mid.seq)[
        , .(proportion, G = G + C)]$G)
    }

    # Calculate G|C content in genome
    win.dt[, paste0("win", win, "_G|C") :=
      genome$count_win_content(chr.names, "G|C", win, regex = TRUE)$`G|C`]
    for (mid.seq in mid.seqs) {
      set(win.dt, j = paste0("win", win, "_G|C_", mid.seq),
          value = genome$count_win_content(chr.names, "G|C", win, mid.seq,
                                           regex = TRUE)$`G|C`)
    }

  }
  fwrite(win.dt, win.file)

  return(win.dt)
}
