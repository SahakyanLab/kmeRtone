removeCaseZone <- function(control.region, case.region, strand.mode,
                           env=parent.frame()) {
  # - Case zone is defined as case region and region on its opposite strand.
  # - The case region is defined as case kmer region plus buffer region.
  # - The case region should be built before using this function.
  # - Strand sensitivity should be resolved i.e. no mixing of * strand and + - strands
  # - Any control regions that overlap WITHIN the damage zone will be removed
  #        e.g. case.zone = 10-20 and control.region = 1-30
  #             the results would be: control.region = 1-10, 20-30
  #     This follows the rule of our control selection where the control is defined to
  #     start AT +80 or -80. The damage zone is -80 --- kmer region --- +80
  #
  # control.region           <string>      A variable name of genomic coordinate table.
  # case.region              <string>      A variable name of genomic coordinate table.
  # env                    <environment>   An environment
  # strand.mode              <string>      "sensitive" or "insensitive". In "sensitive" mode,
  #                                        only "+" or "-" should exist whereas in "insensitive"
  #                                        mode, only "*" should exist.

  if (class(control.region)[1] != "character" | length(control.region) != 1 |
      class(case.region)[1] != "character" | length(case.region) != 1) {
    print()
    stop("Please input genomic.coordinate as a variable name in string format.")
  }
  
  # [1] Build the case zone for strand sensitive mode
  if (strand.mode == "sensitive") {
    ## double the table to reflect case on both strands
    case.zone <- rbind(env[[case.region]][, .(chromosome, start, end, strand)][, strand := "+"],
                       env[[case.region]][, .(chromosome, start, end, strand)][, strand := "-"])
    ## merge
    mergeGenCoordinate("case.zone", environment())
    
  } else if (strand.mode == "insensitive") {
    ## in strand insensitive mode, strand region covers both strands, so it is a case region.
    case.zone <- env[[case.region]][, .(chromosome, start, end, strand)]
    
    ## check if + - exist.
    if (env[[case.region]][strand != "*", .N] > 0 | env[[control.region]][strand != "*", .N] > 0) {
      stop("There is + and - strands in the table. Are you sure this table is in insensitive mode?")
    }
    
  } else stop("Please specify strand.mode either sensitive of insensitive")
  
  # [2] Combine the case zone and control region
  dt <- rbind(env[[control.region]], case.zone)
  setkey(dt, chromosome, strand, start, end)
  
  # [3] Locate continuous region coordinates and assign group
  dt[, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < -1)),
     by = list(chromosome, strand)]
  
  # [4] Make partition of overlapping regions in each group [time consuming]
  # dt <- dt[, .(start = head(unique( sort(c(start, end)) ), -1),
  #              end = unique(sort(c(start, end)))[-1]),
  #          by = list(chromosome, strand, group)]
  
  dt <- dt[, {
    coordinate <- sort( union(start, end) )
    list(start = coordinate[-length(coordinate)], end = coordinate[-1])
  }, by = .(chromosome, strand, group)]
  
  gc()

  # [5] Reorganise the columns
  dt[, group := NULL]
  setcolorder(dt, c("chromosome", "start", "end", "strand"))

  # [6] Remove partitions that are within the case zones
  ## shrink the partitions by one on each of their terminal
  dt[, `:=`(start = start + 1, end = end - 1)]
  
  ## add original un-partitioned case zone to the table
  dt <- rbind(dt, case.zone)
  gc()
  
  ## find overlaps i.e. shrunken partitions within the case zones
  setkey(dt, chromosome, strand, start, end)
  dt[, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < 0)),
     by = list(chromosome, strand)]
  
  ## remove group with more than one partition i.e. the overlaps
  dt <- dt[, if (.N == 1) .(start, end),
            by = .(chromosome, strand, group)]
  gc()
  
  ## expand the partition back
  dt[, `:=`(start = start - 1, end = end + 1)]
  
  # sort the column
  dt[, group := NULL]
  setcolorder(dt, c("chromosome", "start", "end", "strand"))
  
  ########## Alternative to [6] Removing the partitions...
  ### use foverlap function (stable but flagged as experimental by data.table developer)
  ### around the same speed based on human eyes
  # setkey(case.zone, chromosome, strand, start, end)
  # idx <- unique(na.omit(foverlaps(dt, case.zone, type = "within", which = T))$xid)
  # dt1 <- dt[-idx]
  
  # [7] Merge continuous partition
  mergeGenCoordinate("dt", environment())

  env[[control.region]] <- dt
}
