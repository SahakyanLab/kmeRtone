expandGenCoordinate <- function(genomic.coordinate, env=parent.frame(), k) {
  # Expand coordinate to size k. If the coordinate is bigger than k, stop.
  #
  # genomic.coordinate    <string>      A variable name pointing to genomic coordinate <data.table>
  #                                     containing column "chromosome", "start", "end"
  # env                 <environment>   An environment object where the data.table exists.
  
  # Dependency
  #    Package: data.table
  
  if (class(genomic.coordinate)[1] != "character") {
    stop("Please input a variable name in string format.")
  }
  
  # check size
  if (env[[genomic.coordinate]][end-start+1 > k, .N] > 0) {
    stop("The length of coordinates is more than k. The coordinates cannot be expanded.")
  }
  
  # shrink start and end coordinates to point to their midpoint
  env[[genomic.coordinate]][, `:=`(start = start + ((end-start+1)%/%2) - 1,
                            end = end - as.integer((end-start+1)/2 + 0.5) + 1)]
  
  # genomic coordinate length
  len <- env[[genomic.coordinate]][, unique(end - start + 1)]
  
  # if len is not either 1 or 2, something is wrong with the shrinking
  if (!len %in% 1:2 | length(len) > 1) {
    stop("Something is wrong. Contact the developer.")
  }
  
  # scale up to k; if middle point is odd number expand k/2 on both side; if even, expand k/2-1
  env[[genomic.coordinate]][, `:=`(
    start = {
      len <- end-start+1
      start[len == 2] <- start - k/2 + 1
      start[len == 1] <- start - k/2
      start
    },
    end = {
      end[len == 2] <- end + k/2 - 1
      end[len == 1] <- end + k/2
      end
    }
  )]
  
  if (env[[genomic.coordinate]][, unique(end-start+1) > k]) {
    
    message(paste0("\tOdd number of DNA pattern. Size of k is changed from ", k, " to ",
                   genomic.coordinate[, unique(end-start+1)], "."))
    message("k is updated automatically.")
    env$k <- genomic.coordinate[, unique(end-start+1)]
    
  }
}