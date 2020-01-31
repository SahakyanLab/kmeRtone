retrieveCount <- function(genome.name, chromosome, position, width, type, suffix=".csv", FULL.PATH=NULL) {
  # genome        <string> genome used e.g. GRCh37
  # chromosome    <string> e.g. chr1
  # width         <int>
  # type          <string> GC, G, etc.
  # suffix        <string> file suffix e.g. .csv, .csv.gz, etc.
  # FULL.PATH     <string> full path to the reference file. if provided, all other variable
  
  filename <- paste0("data/", type, "/", genome.name, "_", chromosome, "_", type, "_width_", width, suffix)

  first.line <- suppressWarnings(as.numeric(scan(filename, "", nlines = 1, quiet = T)))
  
  # offset position if there is a header in the file
  if (is.na(first.line)) {
    position <- position + 1
  }
  
  position.min <- min(position)
  
  count.range <- scan(filename, integer(), skip = position.min -1, nlines = max(position) - position.min + 1, quiet = T)
  
  # offset the position
  position <- position - position.min + 1
  
  count <- count.range[position]
  
  return(count)
}