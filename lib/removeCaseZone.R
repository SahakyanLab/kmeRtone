removeCaseZone <- function(control.region, case.region, strand.sensitive,
                           env=parent.frame()) {
  # - Case zone is defined as case region and region on its opposite strand.
  # - The case region is defined as case kmer region plus buffer region.
  # - The case region should be built before using this function.
  # - Case zone is case region on both strand.
  # - Here the case zone is referred to as "case.region" after converting - to
  #   + strand.
  # - Strand sensitivity should be resolved i.e. no mixing of * strand and + and
  #   - strands
  # - Any control regions that overlap WITHIN the case zone will be removed
  #        e.g. case.zone = 10-20 and control.region = 1-30
  #             the results would be: control.region = 1-10, 20-30
  #     This follows the rule of our control selection where the control is
  #     defined to start AT +80 or -80.
  #     The case zone is -80 --- kmer region --- +80
  #
  # control.region    <string>  A variable name of genomic coordinate table.
  # case.region       <string>  A variable name of genomic coordinate table.
  # env                <env>    An environment
  # strand.sensitive   <bool>   In strand.sensitive TRUE, only "+" or "-"
  #                             should exist whereas in "insensitive" mode, only
  #                             "*" should exist.
  #
  # Only work on + strand. Control region will be on both strand next to the
  # case zone.

  if (class(control.region)[1] != "character" | length(control.region) != 1 |
      class(case.region)[1] != "character" | length(case.region) != 1) {
    print()
    stop("Please input genomic.coordinate as a variable name in string format.")
  }

  # [1] Combine the case zone and control region
  dt <- rbind(env[[control.region]], env[[case.region]])
  setkey(dt, chromosome, strand, start, end)

  # [2] Locate continuous region coordinates and assign group
  dt[, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < -1)),
     by = list(chromosome, strand)]

  # [3] Make partition of overlapping regions in each group [time consuming]
  # dt <- dt[, .(start = head(unique( sort(c(start, end)) ), -1),
  #              end = unique(sort(c(start, end)))[-1]),
  #          by = list(chromosome, strand, group)]

  dt <- dt[, {
    coordinate <- sort( union(start, end) )
    list(start = coordinate[-length(coordinate)], end = coordinate[-1])
  }, by = .(chromosome, strand, group)]

  gc()

  # [4] Reorganise the columns
  dt[, group := NULL]
  setcolorder(dt, c("chromosome", "start", "end", "strand"))

  # [5] Remove partitions that are within the case zones
  ## shrink the partitions by one on each of their terminal
  dt[, `:=`(start = start + 1, end = end - 1)]

  ## add original un-partitioned case zone to the table
  dt <- rbind(dt, env[[case.region]])
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

  ########## Alternative to [5] Removing the partitions...
  ### use foverlap function (stable but flagged as experimental by data.table
  ### developer)
  ### around the same speed based on human eyes
  # setkey(env[[case.region]], chromosome, strand, start, end)
  # idx <- unique(na.omit(foverlaps(dt, env[[case.region]], type = "within",
  #               which = T))$xid)
  # dt1 <- dt[-idx]

  # [6] Merge continuous partition
  mergeGenCoordinate("dt", environment())

  env[[control.region]] <- dt
}
